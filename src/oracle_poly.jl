#############################################################################
# JuMPeR
# Julia for Mathematical Programming - extension for Robust Optimization
# See http://github.com/IainNZ/JuMPeR.jl
#############################################################################
# PolyhedralOracle
# The default oracle - uses uncertainty bounds and the uncertainty set
# built up with linear constraints to either cut, reformulate, sample or
# some combination of the above.
#############################################################################
using MathProgBase

type PolyhedralOracle <: AbstractOracle
    # The constraints associated with this oracle and the selected mode(s)
    # of operation for each constraint.
    cons::Vector{UncConstraint}
    con_modes::Vector{Dict{Symbol,Bool}}
    con_inds::Dict{Int,Int}
    setup_done::Bool

    # Cutting plane algorithm
    cut_model::Model
    cut_vars::Vector{Variable}
    cut_tol::Float64

    # Reformulation
    num_dualvar::Int
    dual_A::Vector{Vector{(Int,Float64)}}   # Transp. uncertainty matrix
    dual_objs::Vector{Float64}              # Coefficient in dual obj
    dual_vartype::Vector{Symbol}            # >=0, <=0, etc.
    dual_contype::Vector{Symbol}            # Type of new constraint

    # Other options
    debug_printcut::Bool
end
# Default constructor
PolyhedralOracle() = 
    PolyhedralOracle(   UncConstraint[], Dict{Symbol,Bool}[], Dict{Int,Int}(), false,
                        Model(), Variable[], 0.0, # Cutting plane
                        0, Vector{(Int,Float64)}[], Float64[], Symbol[], Symbol[],
                        false)


# registerConstraint
# We must handle this constraint, and the users preferences have been
# communicated through prefs
function registerConstraint(w::PolyhedralOracle, con, ind::Int, prefs)
    con_mode = [:Cut    =>  get(prefs, :prefer_cuts, false), 
                :Reform => !get(prefs, :prefer_cuts, false)]
    push!(w.con_modes, con_mode)
    push!(w.cons, con)
    w.con_inds[ind] = length(w.cons)

    # Extract preferences we care about
    w.debug_printcut = get(prefs, :debug_printcut, false)
    w.cut_tol        = get(prefs, :cut_tol, 1e-6)

    return con_mode
end


# setup
# Now that all modes of operation have been selected, generate the cutting
# plane model and/or the reformulation
function setup(w::PolyhedralOracle, rm::Model)
    w.setup_done && return
    rd = getRobust(rm)
    any_cut = any_reform = false
    for con_mode in w.con_modes
        any_cut    |= con_mode[:Cut]
        any_reform |= con_mode[:Reform]
    end

    # Cutting plane setup
    if any_cut
        # Create an LP that we'll use to solve the cut problem
        # Copy the uncertainty set from the original problem
        w.cut_model.solver   = rd.cutsolver == nothing ? rm.solver : rd.cutsolver
        w.cut_model.numCols  = rd.numUncs
        w.cut_model.colNames = rd.uncNames
        w.cut_model.colLower = rd.uncLower
        w.cut_model.colUpper = rd.uncUpper
        w.cut_model.colCat   = zeros(Int,rd.numUncs)  # TODO: Non-continuous?
        w.cut_vars = [Variable(w.cut_model, i) for i = 1:rd.numUncs]
        for c in rd.uncertaintyset
            newcon = LinearConstraint(AffExpr(), c.lb, c.ub)
            newcon.terms.coeffs = c.terms.coeffs
            newcon.terms.vars   = [w.cut_vars[u.unc] for u in c.terms.vars]
            push!(w.cut_model.linconstr, newcon)
        end
    end

    # Reformulation setup
    if any_reform
        # We have one new variable for every constraint in the uncertainty set
        # We will have one constraint for every uncertainty
        num_dualvar = length(rd.uncertaintyset)
        num_dualcon = rd.numUncs
        dual_A       = [(Int,Float64)[] for i = 1:rd.numUncs]    # Transp. uncertainty matrix
        dual_objs    = zeros(num_dualvar)           # Coefficient in dual obj
        dual_vartype = [:>= for i = 1:num_dualvar]  # >=0, <=0, etc.
        dual_contype = [:>= for i = 1:num_dualcon]  # Type of new constraint

        # Matrix transpose of uncertainty set
        for uncset_i = 1:length(rd.uncertaintyset)
            lhs = rd.uncertaintyset[uncset_i].terms
            for unc_j = 1:length(lhs.vars)
                push!(dual_A[lhs.vars[unc_j].unc], (uncset_i, lhs.coeffs[unc_j]))
            end
            dual_objs[uncset_i] = rhs(rd.uncertaintyset[uncset_i])
            dual_vartype[uncset_i] = sense(rd.uncertaintyset[uncset_i])
        end

        # Add entries for the uncertain bounds
        for unc_i = 1:rd.numUncs
            # If it has a lower bound
            if rd.uncLower[unc_i] != -Inf
                num_dualvar += 1
                push!(dual_A[unc_i], (num_dualvar, 1.0))
                push!(dual_objs, rd.uncLower[unc_i])
                push!(dual_vartype, :>=)
            end
            # If it has an upper bound
            if rd.uncUpper[unc_i] != +Inf
                num_dualvar += 1
                push!(dual_A[unc_i], (num_dualvar, 1.0))
                push!(dual_objs, rd.uncUpper[unc_i])
                push!(dual_vartype, :<=)
            end
            # Now we can handle the sense of the bound in the event
            # one of them is zero. Default is equality.
            if rd.uncLower[unc_i] == 0.0 && rd.uncUpper[unc_i] != 0.0
                dual_contype[unc_i] = :>=
            end
            if rd.uncLower[unc_i] != 0.0 && rd.uncUpper[unc_i] == 0.0
                dual_contype[unc_i] = :<=
            end
            if rd.uncLower[unc_i] != 0.0 && rd.uncUpper[unc_i] != 0.0
                dual_contype[unc_i] = :(==)
            end
        end

        w.num_dualvar   = num_dualvar
        w.dual_A        = dual_A
        w.dual_objs     = dual_objs
        w.dual_vartype  = dual_vartype
        w.dual_contype  = dual_contype 
    end  # end reform preparation
    w.setup_done = true
end


function generateCut(w::PolyhedralOracle, rm::Model, ind::Int, m::Model, cb=nothing, active=false)
    # If not doing cuts for this one, just skip
    con_ind = w.con_inds[ind]
    if !w.con_modes[con_ind][:Cut]
        return 0
    end
    con = w.cons[con_ind]
    rd = getRobust(rm)
    
    master_sol = m.colVal

    # Update the cutting plane problem's objective, and solve
    cut_sense, unc_obj_coeffs, lhs_const = JuMPeR.build_cut_objective(con, master_sol)
    @setObjective(w.cut_model, cut_sense, sum{u[2]*w.cut_vars[u[1]], u=unc_obj_coeffs})
    cut_solve_status = solve(w.cut_model)
    lhs_of_cut = getObjectiveValue(w.cut_model) + lhs_const

    # TEMPORARY: active cut detection
    if active && (
       ((sense(con) == :<=) && (abs(lhs_of_cut - con.ub) <= w.cut_tol)) ||
       ((sense(con) == :>=) && (abs(lhs_of_cut - con.lb) <= w.cut_tol)))
        # Yup its active
        push!(rd.activecuts, w.cut_model.colVal[:])
    end
    
    # Check violation
    if !is_constraint_violated(con, lhs_of_cut, w.cut_tol)
        w.debug_printcut && debug_printcut(rm,m,w,lhs_of_cut,con,nothing)
        return 0  # No violation, no new cut
    end
    
    # Create and add the new constraint
    new_con = JuMPeR.build_certain_constraint(con, w.cut_model.colVal)
    w.debug_printcut && debug_printcut(rm,m,w,lhs_of_cut,con,new_con)
    cb == nothing ? addConstraint(m, new_con) :
                    addLazyConstraint(cb, new_con)
    return 1
end


function generateReform(w::PolyhedralOracle, rm::Model, ind::Int, m::Model)
    # If not doing reform for this one, just skip
    con_ind = w.con_inds[ind]
    if !w.con_modes[con_ind][:Reform]
        return false
    end
    rd = getRobust(rm)

    num_dualvar = w.num_dualvar
    num_dualcon = rd.numUncs
    dual_A      = w.dual_A
    dual_objs   = w.dual_objs
    dual_vartype= w.dual_vartype
    dual_contype= w.dual_contype 
    dual_rhs    = [AffExpr() for i = 1:num_dualcon]
    new_lhs     = AffExpr()  # This will be the new constraint
    new_rhs     = rhs(w.cons[con_ind])
    orig_lhs    = w.cons[con_ind].terms

    # Collect certain terms
    for term_i = 1:length(orig_lhs.vars)
        if orig_lhs.coeffs[term_i].constant != 0
            # Constant part
            push!(new_lhs, orig_lhs.coeffs[term_i].constant,
                          Variable(m, orig_lhs.vars[term_i].col))
        end
    end
    
    # Collect coefficients for each uncertainty present
    for term_i = 1:length(orig_lhs.coeffs)
        for unc_j = 1:length(orig_lhs.coeffs[term_i].coeffs)
            push!(dual_rhs[orig_lhs.coeffs[term_i].vars[unc_j].unc],
                    orig_lhs.coeffs[term_i].coeffs[unc_j],
                    Variable(m, orig_lhs.vars[term_i].col))
        end
    end
        for unc_j = 1:length(orig_lhs.constant.coeffs)
            dual_rhs[orig_lhs.constant.vars[unc_j].unc] +=
                                orig_lhs.constant.coeffs[unc_j]
        end

    # Create dual variables for this contraint, and add to new LHS
    dual_vars = Variable[]
    for dual_i = 1:num_dualvar
        # Constraint is less-than, so "maximize"
        var_name = "_µ$(con_ind)_$(dual_i)"
        if sense(w.cons[con_ind]) == :<=
            if dual_vartype[dual_i] == :<=      # LE  ->  >= 0
                push!(dual_vars, Variable(m,0,+Inf,0,var_name))
            elseif dual_vartype[dual_i] == :>=  # GE  ->  <= 0
                push!(dual_vars, Variable(m,-Inf,0,0,var_name))
            elseif dual_vartype[dual_i] == :(==)  # EQ  ->  free
                push!(dual_vars, Variable(m,-Inf,+Inf,0,var_name))
            end
        end
        # Constraint is gerater-than, so "minimize"
        if sense(w.cons[con_ind]) == :>=
            if dual_vartype[dual_i] == :>=      # GE  ->  >= 0
                push!(dual_vars, Variable(m,0,+Inf,0,var_name))
            elseif dual_vartype[dual_i] == :<=  # LE  ->  <= 0
                push!(dual_vars, Variable(m,-Inf,0,0,var_name))
            elseif dual_vartype[dual_i] == :(==)  # EQ  ->  free
                push!(dual_vars, Variable(m,-Inf,+Inf,0,var_name))
            end
        end

        push!(new_lhs, dual_objs[dual_i], dual_vars[dual_i])
    end

    # Add the new constraint which replaces the original constraint
    sense(w.cons[con_ind]) == :(<=) && addConstraint(m, new_lhs <= new_rhs)
    sense(w.cons[con_ind]) == :(>=) && addConstraint(m, new_lhs >= new_rhs)

    # Add the additional new constraints
    for unc_i = 1:rd.numUncs
        new_lhs = AffExpr()
        for pair in dual_A[unc_i]
            push!(new_lhs, pair[2], dual_vars[pair[1]])
        end
        ucontype = sense(w.cons[con_ind])
        dualtype = dual_contype[unc_i]
        if      dualtype == :(==)
            addConstraint(m, new_lhs == dual_rhs[unc_i])
        elseif (dualtype == :(<=) && ucontype == :(<=)) ||
               (dualtype == :(>=) && ucontype == :(>=))
            addConstraint(m, new_lhs <= dual_rhs[unc_i])
        else
            addConstraint(m, new_lhs >= dual_rhs[unc_i])
        end
    end

    return true
end


function debug_printcut(rm,m,w,lhs,con,new_con)
    println("BEGIN DEBUG :debug_printcut")
    convert_model!(con, rm)
    println("  Constraint:  ", con)
    convert_model!(con, m)
    println("  Master sol:  ")
    for j in 1:length(m.colNames)
        println("    ", m.colNames[j], "  ", m.colVal[j])
    end
    print(w.cut_model)
    #println("  Solve status ", cut_solve_status)
    println("  Cut sol:     ")
    for j in 1:length(w.cut_model.colNames)
        println("    ", w.cut_model.colNames[j], "  ", w.cut_model.colVal[j])
    end
    println("  OrigLHS val: ", lhs)
    println("  Sense:       ", sense(con))
    println("  con.lb/ub:   ", con.lb, "  ", con.ub)
    println("  new con:  ", new_con)
    println("END DEBUG   :debug_printcut")
end