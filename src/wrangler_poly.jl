#############################################################################
# JuMP
# An algebraic modelling langauge for Julia
# See http://github.com/JuliaOpt/JuMP.jl
#############################################################################
# PolyhedralWrangler
# The default wrangler - uses uncertainty bounds and the uncertainty set
# built up with linear constraints to either cut, reformulate, sample or
# some combination of the above.
#############################################################################

type PolyhedralWrangler <: AbstractWrangler
    # The constraints associated with this wrangler and the selected mode(s)
    # of operation for each constraint.
    cons::Vector{UncConstraint}
    con_modes::Vector{Dict{Symbol,Bool}}
    con_inds::Dict{Int,Int}
    any_reform::Bool
    any_cut::Bool
    setup_done::Bool

    # Cutting plane algorithm
    cut_model::Model
    cut_vars::Vector{Variable}
    var_maps    # Mappings from master variable values to cutting objective
    const_maps  # The constant parts (not affected by uncertainty)

    # Reformulation
    num_dualvar::Int
    dual_A::Vector{Vector{(Int,Float64)}}   # Transp. uncertainty matrix
    dual_objs::Vector{Float64}              # Coefficient in dual obj
    dual_vartype::Vector{Symbol}            # >=0, <=0, etc.
    dual_contype::Vector{Symbol}            # Type of new constraint
end
# Default constructor
PolyhedralWrangler() = 
    PolyhedralWrangler(UncConstraint[], Dict{Symbol,Bool}[], Dict{Int,Int}(),
                       false, false, false,
                       Model(), Variable[], Any[], Any[],  # Cutting plane
                       0, Vector{(Int,Float64)}[], Float64[], Symbol[], Symbol[])


# registerConstraint
# We must handle this constraint, and the users preferences have been
# communicated through prefs
function registerConstraint(w::PolyhedralWrangler, con, ind::Int, prefs)
    push!(w.cons, con)
    w.con_inds[ind] = length(w.cons)
    con_mode = Dict{Symbol,Bool}()
    if :prefer_cuts in keys(prefs) && prefs[:prefer_cuts]
        con_mode = [:Cut => true, :Reform => false,
                    :Sample => (:samples in keys(prefs)) ? 
                                prefs[:samples] > 0 : false  ]
        w.any_cut = true
    else  # Default to reformulation
        con_mode = [:Cut => false, :Reform => true, :Sample => false]
        w.any_reform = true
    end
    push!(w.con_modes, con_mode)
    return con_mode
end


# setup
# Now that all modes of operation have been selected, generate the cutting
# plane model and/or the reformulation
function setup(w::PolyhedralWrangler, rm::Model)
    if w.setup_done
        return
    end
    rd = getRobust(rm)

    # Cutting plane setup
    if w.any_cut
        # Create an LP that we'll use to solve the cut problem
        # Copy the uncertainty set from the original problem
        w.cut_model.solver   = rm.solver
        w.cut_model.numCols  = rd.numUncs
        w.cut_model.colNames = rd.uncNames
        w.cut_model.colLower = rd.uncLower
        w.cut_model.colUpper = rd.uncUpper
        w.cut_model.colCat   = zeros(rd.numUncs)  # TODO: Non-continuous?
        w.cut_vars = [Variable(w.cut_model, i) for i = 1:rd.numUncs]
        for c in rd.uncertaintyset
            newcon = LinearConstraint(AffExpr(), c.lb, c.ub)
            newcon.terms.coeffs = c.terms.coeffs
            newcon.terms.vars   = [w.cut_vars[u.unc] for u in c.terms.uncs]
            push!(w.cut_model.linconstr, newcon)
        end

        # For each constraint we are doing cuts for, build up a mapping
        # from the master solution to the cut objective
        # w.var_map: key = column numbers, value = [(unc_ind, coeff)]
        # w.constant_coeffs = [(unc_ind, coeff)]
        for con_ind = 1:length(w.cons)
            # Don't build a map if not doing cuts
            if !w.con_modes[con_ind][:Cut]
                push!(w.var_maps, nothing)
                push!(w.const_maps, nothing)
                continue
            end
            
            con_lhs = w.cons[con_ind].terms
            num_var = length(con_lhs.coeffs)  # Number of variables in constraint
            var_map = Dict{Int,Vector{(Int,Float64)}}()
            for i = 1:num_var                               # For every (unc_aff, var) pair
                col = con_lhs.vars[i].col                   # Extract column number [key]
                var_map[col] = (Int,Float64)[]              # Init array of tuples
                num_unc = length(con_lhs.coeffs[i].uncs)    # Num unc in unc_aff
                sizehint(var_map[col], num_unc)
                for j = 1:num_unc                           # For every unc in unc_aff
                    unc   = con_lhs.coeffs[i].uncs[j].unc   # The index of this uncertain
                    coeff = con_lhs.coeffs[i].coeffs[j]     # The coeff on this uncertain
                    push!(var_map[col], (unc,coeff))        # Add to mapping
                end
            end
            push!(w.var_maps, var_map)
            
            # Including the constant term of the constraint
            const_map = (Int,Float64)[]
                for j = 1:length(con_lhs.constant.uncs)
                    unc   = con_lhs.constant.uncs[j].unc
                    coeff = con_lhs.constant.coeffs[j]
                    push!(const_map, (unc,coeff))
                end
            push!(w.const_maps, const_map)
        end  # Next constraint
    end  # end cut preparation

    # Reformulation setup
    if w.any_reform
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
            for unc_j = 1:length(lhs.uncs)
                push!(dual_A[lhs.uncs[unc_j].unc], (uncset_i, lhs.coeffs[unc_j]))
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


function generateCut(w::PolyhedralWrangler, rm::Model, ind::Int, m::Model)
    
    # If not doing cuts for this one, just skip
    con_ind = w.con_inds[ind]
    if !w.con_modes[con_ind][:Cut]
        return 0
    end
    rd = getRobust(rm)

    master_sol = m.colVal

    # Objective of cut problem depends on sense of constraint
    con = w.cons[con_ind]
    w.cut_model.objSense = sense(con) == :<= ? :Max : :Min

    # Shove the master solution into the objective using our map
    num_uncs = w.cut_model.numCols
    unc_coeffs = zeros(num_uncs)
    for key in keys(w.var_maps[con_ind])        # For every var in con
        for pair in w.var_maps[con_ind][key]    # For every unc in front of that var
            unc_coeffs[pair[1]] += pair[2]*master_sol[key]
        end
    end
        for pair in w.const_maps[con_ind]
            unc_coeffs[pair[1]] += pair[2]
        end
    @setObjective(w.cut_model, w.cut_model.objSense, sum{unc_coeffs[i]*w.cut_vars[i], i = 1:num_uncs})

    # Solve cutting plane problem
    solve(w.cut_model)

    # Calculate cut LHS
    lhs = 0.0
    orig_lhs = con.terms
    # Variable part
    num_vars = length(orig_lhs.vars)
    for var_ind = 1:num_vars
        coeff   = orig_lhs.coeffs[var_ind]                  # The unc expr in front of var
        col_val = master_sol[orig_lhs.vars[var_ind].col]    # The value of this x
        lhs    += coeff.constant * col_val                 # The constant part for this x
        for unc_ind = 1:length(coeff.uncs)
            coeff_unc   = coeff.uncs[unc_ind]
            coeff_coeff = coeff.coeffs[unc_ind]
            lhs += w.cut_model.colVal[coeff_unc.unc] * coeff_coeff[unc_ind] * col_val
        end
    end
    # Non variable part
    coeff = orig_lhs.constant
    lhs  += coeff.constant
        for unc_ind = 1:length(coeff.uncs)
            coeff_unc   = coeff.uncs[unc_ind]
            coeff_coeff = coeff.coeffs[unc_ind]
            lhs += w.cut_model.colVal[coeff_unc.unc] * coeff_coeff[unc_ind]
        end

    # Check violation
    if ((sense(con) == :<=) && (lhs <= con.ub + 1e-6)) ||
       ((sense(con) == :>=) && (lhs >= con.lb - 1e-6))
        return 0  # No violation, no new cut
    end
    
    # Now add that solution back in
    # TODO: Build map for this too?
    new_lhs = AffExpr(orig_lhs.vars,
                      [orig_lhs.coeffs[i].constant for i in 1:num_vars],
                      orig_lhs.constant.constant)
    # Variable part
    for var_ind = 1:num_vars
        coeff   = orig_lhs.coeffs[var_ind]
        for unc_ind = 1:length(coeff.uncs)
            coeff_unc   = coeff.uncs[unc_ind]
            coeff_coeff = coeff.coeffs[unc_ind]
            new_lhs.coeffs[var_ind] += w.cut_model.colVal[coeff_unc.unc] * coeff_coeff[unc_ind]
        end
    end
    # Non variable part
        coeff = orig_lhs.constant
        for unc_ind = 1:length(coeff.uncs)
            coeff_unc = coeff.uncs[unc_ind]
            coeff_coeff = coeff.coeffs[unc_ind]
            new_lhs.constant += w.cut_model.colVal[coeff_unc.unc] * coeff_coeff[unc_ind]
        end

    if sense(con) == :<=
        @addConstraint(m, new_lhs <= con.ub)
    else
        @addConstraint(m, new_lhs >= con.lb)
    end
    return 1
end


function generateReform(w::PolyhedralWrangler, rm::Model, ind::Int, m::Model)
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
            push!(new_lhs.vars, Variable(m, orig_lhs.vars[term_i].col))
            push!(new_lhs.coeffs, orig_lhs.coeffs[term_i].constant)
        end
    end
    
    # Collect coefficients for each uncertainty present
    for term_i = 1:length(orig_lhs.coeffs)
        for unc_j = 1:length(orig_lhs.coeffs[term_i].coeffs)
            dual_rhs[orig_lhs.coeffs[term_i].uncs[unc_j].unc] +=
                                Variable(m, orig_lhs.vars[term_i].col) *
                                    orig_lhs.coeffs[term_i].coeffs[unc_j]
        end
    end
        for unc_j = 1:length(orig_lhs.constant.coeffs)
            dual_rhs[orig_lhs.constant.uncs[unc_j].unc] +=
                                orig_lhs.constant.coeffs[unc_j]
        end

    # Create dual variables for this contraint, and add to new LHS
    dual_vars = Variable[]
    for dual_i = 1:num_dualvar
        # Constraint is less-than, so "maximize"
        if sense(w.cons[con_ind]) == :<=
            if dual_vartype[dual_i] == :<=      # LE  ->  >= 0
                push!(dual_vars, Variable(m,0,+Inf,0))
            elseif dual_vartype[dual_i] == :>=  # GE  ->  <= 0
                push!(dual_vars, Variable(m,-Inf,0,0))
            elseif dual_vartype[dual_i] == :(==)  # EQ  ->  free
                push!(dual_vars, Variable(m,-Inf,+Inf,0))
            end
        end
        # Constraint is gerater-than, so "minimize"
        if sense(w.cons[con_ind]) == :>=
            if dual_vartype[dual_i] == :>=      # GE  ->  >= 0
                push!(dual_vars, Variable(m,0,+Inf,0))
            elseif dual_vartype[dual_i] == :<=  # LE  ->  <= 0
                push!(dual_vars, Variable(m,-Inf,0,0))
            elseif dual_vartype[dual_i] == :(==)  # EQ  ->  free
                push!(dual_vars, Variable(m,-Inf,+Inf,0))
            end
        end

        new_lhs += dual_objs[dual_i]*dual_vars[dual_i]
    end

    # Add the new constraint which replaces the original constraint
    if sense(w.cons[con_ind]) == :<=
        addConstraint(m, new_lhs <= new_rhs)
    elseif sense(w.cons[con_ind]) == :>=
        addConstraint(m, new_lhs >= new_rhs)
    end

    # Add the additional new constraints
    for unc_i = 1:rd.numUncs
        new_lhs = AffExpr()
        for pair in dual_A[unc_i]
            new_lhs += pair[2] * dual_vars[pair[1]]
        end
        if dual_contype[unc_i] == :(==)
            addConstraint(m, new_lhs == dual_rhs[unc_i])
        elseif dual_contype[unc_i] == :<=
            addConstraint(m, new_lhs <= dual_rhs[unc_i])
        elseif dual_contype[unc_i] == :>=
            addConstraint(m, new_lhs >= dual_rhs[unc_i])
        end
    end

    return true
end

