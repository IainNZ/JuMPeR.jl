#############################################################################
# JuMPeR
# Julia for Mathematical Programming - extension for Robust Optimization
# See http://github.com/IainNZ/JuMPeR.jl
#############################################################################
# GeneralOracle
# The default oracle - uses uncertainty bounds and the uncertainty set
# built up with linear constraints and ellipses and either use cutting planes
# or reformulate to obtain a deterministic problem.
#############################################################################

type GeneralOracle <: AbstractOracle
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

    # Reformulation (see setup for comments)
    num_dualvar::Int
    dual_A::Vector{Vector{(Int,Float64)}}
    dual_objs::Vector{Float64}
    dual_vartype::Vector{Symbol}
    dual_contype::Vector{Symbol}
    dual_ell_rhs_idxs::Vector{Vector{Int}}
    dual_ell_lhs_idxs::Vector{Int}

    # Other options
    debug_printcut::Bool
    debug_printreform::Bool
end
# Default constructor
GeneralOracle() = 
    GeneralOracle(  UncConstraint[], Dict{Symbol,Bool}[], Dict{Int,Int}(), false,
                    Model(), Variable[], 0.0, # Cutting plane
                    0, Vector{(Int,Float64)}[], Float64[], Symbol[], Symbol[], 
                    Vector{Int}[], Int[],
                    false, false)


# registerConstraint
# We must handle this constraint, and the users preferences have been
# communicated through prefs
function registerConstraint(w::GeneralOracle, con, ind::Int, prefs)
    con_mode = [:Cut    =>  get(prefs, :prefer_cuts, false), 
                :Reform => !get(prefs, :prefer_cuts, false)]
    push!(w.con_modes, con_mode)
    push!(w.cons, con)
    w.con_inds[ind] = length(w.cons)

    # Extract preferences we care about
    w.debug_printcut    = get(prefs, :debug_printcut, false)
    w.debug_printreform = get(prefs, :debug_printreform, false)
    w.cut_tol           = get(prefs, :cut_tol, 1e-6)

    return con_mode
end


# setup
# Now that all modes of operation have been selected, generate the cutting
# plane model and/or the reformulation
function setup(w::GeneralOracle, rm::Model)
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
        w.cut_model.colCat   = rd.uncCat
        w.cut_vars = [Variable(w.cut_model, i) for i = 1:rd.numUncs]
        # Polyhedral constraints
        for c in rd.uncertaintyset
            newcon = LinearConstraint(AffExpr(), c.lb, c.ub)
            newcon.terms.coeffs = c.terms.coeffs
            newcon.terms.vars   = [w.cut_vars[u.unc] for u in c.terms.vars]
            push!(w.cut_model.linconstr, newcon)
        end
        # Ellipse constraints
        # Take || Fu + g || <= Gamma and rewrite as
        # y = Fu+g, y^T y <= Gamma^2
        for el_c in rd.normconstraints
            yty = QuadExpr()
            num_terms, num_uncs = size(el_c.F)
            for term_ind = 1:num_terms
                # Create new variable y_i
                y   = Variable(w.cut_model,-Inf,Inf,JuMP.CONTINUOUS,"_el_$term_ind")
                Fug = AffExpr(Variable[w.cut_vars[i] for i in el_c.u],
                              Float64[el_c.F[term_ind,i] for i in 1:num_uncs],
                              el_c.g[term_ind])
                addConstraint(w.cut_model, y == Fug)
                push!(yty.qvars1, y)
                push!(yty.qvars2, y)
                push!(yty.qcoeffs, 1.0)
            end
            addConstraint(w.cut_model, yty <= el_c.Gamma^2)
        end
    end

    # Reformulation setup
    if any_reform
        #####################################################################
        # We have 
        # - one new variable for every constraint in the uncertainty set
        # - one new variable for every term in each ellipse
        # - one new variable for each ellipse in general
        # We have one constraint for every uncertainty.
        #
        # Statement of primal-dual pair for the one ellipse case:
        # PRIMAL
        # max  cx.x  + cy.y
        #  st  Ax.x  + Ay.y   ==  b
        #          ||Fy+g||2  <=  Gamma
        #
        # DUAL
        # min   b.pi +  g.bz + Gamma*bt
        #  st Ax'.pi          == cx
        #     Ay'.pi - F'.bz  == cy
        #         pi free
        #             (bt,bz) in Q^n+1
        #
        # In this setup phase we build the structure of the dual but do not
        # attach it to the model. Rather, for each constraint we will spawn
        # a new set of variables and constraints according to the structure
        # we determine here.
        #####################################################################
        num_dualvar = length(rd.uncertaintyset)
        for el_c in rd.normconstraints
            num_dualvar += length(el_c.u) + 1
        end
        num_dualcon = rd.numUncs
        # Sparse representation of the transpose of the matrix defining
        # the linear constraints
        dual_A       = [(Int,Float64)[] for i = 1:rd.numUncs]
        # Objective coefficient in dual obj
        dual_objs    = zeros(num_dualvar)
        # Primal constr. sense -> dual variable sense, >=0, <=0, free (:==)
        # Initially it is just the primal constraint sense, but might need
        # to be flipped later
        dual_vartype = [:>= for i = 1:num_dualvar]
        # Dual constr. sense
        dual_contype = [:>= for i = 1:num_dualcon]

        # First, use linear constraints to
        # - form matrix transpose of linear constraints
        # - set objective coefficients for linear constraint duals
        # - set variable types for linear constraint duals
        for uncset_i = 1:length(rd.uncertaintyset)
            lhs = rd.uncertaintyset[uncset_i].terms
            for unc_j = 1:length(lhs.vars)
                push!(dual_A[lhs.vars[unc_j].unc], (uncset_i, lhs.coeffs[unc_j]))
            end
            dual_objs[uncset_i]    = rhs(rd.uncertaintyset[uncset_i])
            dual_vartype[uncset_i] = sense(rd.uncertaintyset[uncset_i])
        end

        # Next, use the ellipses to set objective coefficients of the duals
        # relating to the ellipses
        start_index = length(rd.uncertaintyset)
        for el_c in rd.normconstraints
            # Duals corresponding to terms in ellipse
            dual_ell_rhs_idx = Int[]
            for term_ind = 1:length(el_c.g)
                # The index of this dual variable
                dual_idx = start_index+term_ind
                # Objective coefficient is constant part of this
                # term in the ellipse
                dual_objs[dual_idx]     = el_c.g[term_ind]
                dual_vartype[dual_idx]  = :Qrhs
                # Store the indices for this ellipse so we inject them
                # into constraints later
                push!(dual_ell_rhs_idx, dual_idx)
            end
            push!(w.dual_ell_rhs_idxs, dual_ell_rhs_idx)

            # Dual corresponding to the RHS dual
            lhs_idx = start_index+length(el_c.g)+1
            dual_objs[lhs_idx]      = el_c.Gamma
            dual_vartype[lhs_idx]   = :Qlhs
            push!(w.dual_ell_lhs_idxs, lhs_idx)
            start_index += length(el_c.g) + 1
        end

        # We will handle bounds on the coefficients as linear constraints,
        # which means adding to the matrix transpose and objective
        for unc_i = 1:rd.numUncs
            # If it has a lower bound...
            if rd.uncLower[unc_i] != -Inf
                num_dualvar += 1
                push!(dual_A[unc_i], (num_dualvar, 1.0))
                push!(dual_objs,     rd.uncLower[unc_i])
                push!(dual_vartype,  :(>=))
            end
            # If it has an upper bound
            if rd.uncUpper[unc_i] != +Inf
                num_dualvar += 1
                push!(dual_A[unc_i], (num_dualvar, 1.0))
                push!(dual_objs,     rd.uncUpper[unc_i])
                push!(dual_vartype,  :(<=))
            end
            # Now we just treat the variable as being free
            dual_contype[unc_i] = :(==)
        end

        w.num_dualvar   = num_dualvar
        w.dual_A        = dual_A
        w.dual_objs     = dual_objs
        w.dual_vartype  = dual_vartype
        w.dual_contype  = dual_contype

        if w.debug_printreform
            println("BEGIN DEBUG :debug_printreform")
            println("Num dual var ", num_dualvar)
            println("Objective")
            dump(dual_objs)
            println("dual_A")
            println(dual_A)
            println("lhs,rhs idx")
            dump(w.dual_ell_lhs_idxs)
            dump(w.dual_ell_rhs_idxs)
            println("END DEBUG   :debug_printreform")
        end
    end  # end reformulation preparation
    w.setup_done = true
end


function generateCut(w::GeneralOracle, rm::Model, ind::Int, m::Model, cb=nothing, active=false)
    # If not doing cuts for this one, just skip
    con_ind = w.con_inds[ind]
    if !w.con_modes[con_ind][:Cut] && !active
        return 0
    end
    con = w.cons[con_ind]
    rd = getRobust(rm)
    
    master_sol = m.colVal

    # Update the cutting plane problem's objective, and solve
    cut_sense, unc_obj_coeffs, lhs_const = JuMPeR.build_cut_objective(con, master_sol)
    @setObjective(w.cut_model, cut_sense, sum{u[2]*w.cut_vars[u[1]], u=unc_obj_coeffs})
    cut_solve_status = solve(w.cut_model, suppress_warnings=true)
    cut_solve_status != :Optimal && error("Cutting plane problem infeasible or unbounded!")
    lhs_of_cut = getObjectiveValue(w.cut_model) + lhs_const

    # TEMPORARY: active cut detection
    if active
        push!(rd.activecuts[ind], 
            cut_to_scen(w.cut_model.colVal, 
                check_cut_status(con, lhs_of_cut, w.cut_tol) == :Active))
    end
    
    # Check violation
    if check_cut_status(con, lhs_of_cut, w.cut_tol) != :Violate
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


function generateReform(w::GeneralOracle, rm::Model, ind::Int, master::Model)
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
    # These are the objective coefficients of the cutting plane problem
    # that are the right-hand-side of the dual problem.
    dual_rhs    = [AffExpr() for i = 1:num_dualcon]
    # This is the reformulated "main" constraint that will replace the
    # existing uncertain constraint
    new_lhs     = AffExpr()

    # We will do all reformulation as if the constraint is a <=
    # constraint. This necessitates mulitplying through by -1
    # if the original constraint is not of this form
    sign_flip   = sense(w.cons[con_ind]) == :(<=) ? +1.0 : -1.0
    # The right-hand-side of the new constraint
    new_rhs     = rhs(w.cons[con_ind]) * sign_flip
    # The uncertain constraint left-hand-side
    orig_lhs    = w.cons[con_ind].terms

    # First collect the certain terms of the uncertain constraint - they
    # won't take part in the reformulation, so we can append them to the
    # new LHS directly
    num_lhs_terms = length(orig_lhs.vars)
    for term_i = 1:num_lhs_terms
        var_col = orig_lhs.vars[term_i].col
        if orig_lhs.coeffs[term_i].constant != 0.0
            push!(new_lhs,  orig_lhs.coeffs[term_i].constant * sign_flip,
                            Variable(master, var_col) )
        end
    end
    
    # We have terms (a^T u + b) x_i, we now need to get (c^T x) u_j
    # The "c^T x" will form the new right-hand-side
    # At same time, we'll check for integral uncertainties, which
    # are a not allowed in the reformulation
    for term_i = 1:num_lhs_terms
        var_col     = orig_lhs.vars[term_i].col
        term_coeff  = orig_lhs.coeffs[term_i]
        for coeff_term_j = 1:length(term_coeff.coeffs)
            coeff = term_coeff.coeffs[coeff_term_j]
            unc   = term_coeff.vars[coeff_term_j].unc
            rd.uncCat[unc] != JuMP.CONTINUOUS && 
                error("No integer uncertainties allowed in reformulation!")
            push!(dual_rhs[unc], coeff * sign_flip, 
                                 Variable(master, var_col))
        end
    end
    # We also need the standalone (a^T u) not related to any variable
        for const_term_j = 1:length(orig_lhs.constant.coeffs)
            coeff = orig_lhs.constant.coeffs[const_term_j]
            unc   = orig_lhs.constant.vars[const_term_j].unc
            rd.uncCat[unc] != JuMP.CONTINUOUS && 
                error("No integer uncertainties allowed in reformulation!")
            dual_rhs[unc].constant += coeff * sign_flip
        end

    # The right-hand-side is the only thing unique to this particular
    # constraint. We now need to create dual variables for this specific
    # contraint and add them to the new constraint's LHS
    dual_vars = Variable[]
    for ind = 1:num_dualvar
        # Determine bounds for case where contraint is <= (maximize)
        # then flip if otherwise
        var_name = "_µ$(con_ind)_$(ind)"
        lbound = -Inf
        ubound = +Inf
        vt = dual_vartype[ind]
        # Less-than and LHS of cone -->  v >= 0
        if vt == :(<=) || vt == :Qlhs
            lbound = 0
        end
        # Greater-than -->  v <= 0
        if vt == :(>=)
            ubound = 0
        end
        # Equality and RHS of cone -->  free
        new_v = Variable(master,lbound,ubound,JuMP.CONTINUOUS,var_name)
        push!(dual_vars, new_v)
        push!(new_lhs, dual_objs[ind], new_v)
    end

    w.debug_printreform && println("DEBUG new_lhs ", new_lhs)

    # Add the new constraint which "replaces" the original constraint
    addConstraint(master, new_lhs <= new_rhs)

    # Add the additional new constraints
    # - If the uncertainty wasn't involved in an ellipse, its just dual_A
    # - If it was, we need to add in additional terms (from .F)
    for unc = 1:rd.numUncs
        new_lhs = AffExpr()
        # Add the linear constraint segment
        for pair in dual_A[unc]
            push!(new_lhs, pair[2], dual_vars[pair[1]])
        end

        # Check if this uncertain is in an ellipse
        ell_idx = 0
        for el_c in rd.normconstraints
            ell_idx += 1
            el_matrix_col = 0
            for el_term = 1:length(el_c.u)
                if el_c.u[el_term] == unc
                    el_matrix_col = el_term
                end
            end
            el_matrix_col == 0 && continue
            F_slice = el_c.F[:,el_matrix_col]
            for el_matrix_row = 1:length(F_slice)
                push!(new_lhs, -F_slice[el_matrix_row], 
                      dual_vars[w.dual_ell_rhs_idxs[ell_idx][el_matrix_row]] )
            end
        end

        dualtype = dual_contype[unc]
        if     dualtype == :(==)
            addConstraint(master, new_lhs == dual_rhs[unc])
        elseif dualtype == :(<=)
            addConstraint(master, new_lhs <= dual_rhs[unc])
        else
            addConstraint(master, new_lhs >= dual_rhs[unc])
        end
    end

    # Finally we impose that the duals relating to the ellipse
    # live in a cone
    ell_idx = 0
    for el_c in rd.normconstraints
        ell_idx += 1
        beta_t = dual_vars[w.dual_ell_lhs_idxs[ell_idx]]
        beta_zs = dual_vars[w.dual_ell_rhs_idxs[ell_idx]]
        addConstraint(master, beta_t^2 >= sum([beta_z^2 for beta_z in beta_zs]))
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
