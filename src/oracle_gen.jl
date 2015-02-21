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
    use_cuts::Bool

    # Cutting plane algorithm
    cut_model::Model
    cut_vars::Vector{Variable}
    cut_tol::Float64

    # Reformulation (see setup() for details)
    num_dualvar::Int
    dual_A::Vector{Vector{(Int,Float64)}}
    dual_objs::Vector{Float64}
    dual_vartype::Vector{Symbol}
    dual_contype::Vector{Symbol}
    dual_ell_rhs_idxs::Vector{Vector{Int}}
    dual_ell_lhs_idxs::Vector{Int}

    # Options
    debug_printcut::Bool
    debug_printreform::Bool
end
# Default constructor
GeneralOracle() = 
    GeneralOracle(  false,
                    Model(), Variable[], 0.0, # Cutting plane
                    0, Vector{(Int,Float64)}[], Float64[], Symbol[], Symbol[], 
                    Vector{Int}[], Int[],
                    false, false)


# registerConstraint
# We must handle this constraint, and the users preferences have been
# communicated through prefs. We don't need to take any action here.
registerConstraint(gen::GeneralOracle, rm::Model, ind::Int, prefs) = nothing


# setup
# Generate the cutting plane model or precompute the reformulation structure.
function setup(gen::GeneralOracle, rm::Model, prefs)

    # Extract preferences we care about
    gen.use_cuts          = get(prefs, :prefer_cuts, false)
    gen.cut_tol           = get(prefs, :cut_tol, 1e-6)
    gen.debug_printcut    = get(prefs, :debug_printcut, false)
    gen.debug_printreform = get(prefs, :debug_printreform, false)

    rd = getRobust(rm)

    # Cutting plane setup
    if gen.use_cuts || prefs[:active_cuts]
        # Create an LP that we'll use to solve the cut problem
        # Copy the uncertainty set from the original problem
        gen.cut_model = Model()
        gen.cut_model.solver   = isa(rd.cutsolver,JuMP.UnsetSolver) ? rm.solver : rd.cutsolver
        gen.cut_model.numCols  = rd.numUncs
        gen.cut_model.colNames = rd.uncNames
        gen.cut_model.colLower = rd.uncLower
        gen.cut_model.colUpper = rd.uncUpper
        gen.cut_model.colCat   = rd.uncCat
        gen.cut_vars = [Variable(gen.cut_model, i) for i = 1:rd.numUncs]
        # Polyhedral constraints
        for c in rd.uncertaintyset
            newcon = LinearConstraint(AffExpr(), c.lb, c.ub)
            newcon.terms.coeffs = c.terms.coeffs
            newcon.terms.vars   = [gen.cut_vars[u.unc] for u in c.terms.vars]
            push!(gen.cut_model.linconstr, newcon)
        end
        # Ellipse constraints
        # Take || Fu + g || <= Gamma and rewrite as
        # y = Fu+g, y^T y <= Gamma^2
        for el_c in rd.normconstraints
            yty = QuadExpr()
            num_terms, num_uncs = size(el_c.F)
            for term_ind = 1:num_terms
                # Create new variable y_i
                y   = Variable(gen.cut_model,-Inf,Inf,:Cont,"_el_$term_ind")
                Fug = AffExpr(Variable[gen.cut_vars[i] for i in el_c.u],
                              Float64[el_c.F[term_ind,i] for i in 1:num_uncs],
                              el_c.g[term_ind])
                addConstraint(gen.cut_model, y == Fug)
                push!(yty.qvars1, y)
                push!(yty.qvars2, y)
                push!(yty.qcoeffs, 1.0)
            end
            # Normally we'd do
            # addConstraint(gen.cut_model, yty <= el_c.Gamma^2)
            # But unfortunately this isn't a SOC constraint so ECOS
            # won't like it (Gurobi can handle it). We just need
            # an auxiliary variable set to equal the RHS.
            rhs_aux = Variable(gen.cut_model,0,Inf,:Cont,"_el_rhs_aux")
            @addConstraint(gen.cut_model, rhs_aux == el_c.Gamma)
            @addConstraint(gen.cut_model, yty <= rhs_aux^2)
        end
    end

    if !gen.use_cuts
        #####################################################################
        # Reformulation setup
        #
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
            push!(gen.dual_ell_rhs_idxs, dual_ell_rhs_idx)

            # Dual corresponding to the RHS dual
            lhs_idx = start_index+length(el_c.g)+1
            dual_objs[lhs_idx]      = el_c.Gamma
            dual_vartype[lhs_idx]   = :Qlhs
            push!(gen.dual_ell_lhs_idxs, lhs_idx)
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

        gen.num_dualvar   = num_dualvar
        gen.dual_A        = dual_A
        gen.dual_objs     = dual_objs
        gen.dual_vartype  = dual_vartype
        gen.dual_contype  = dual_contype

        if gen.debug_printreform
            println("BEGIN DEBUG :debug_printreform")
            println("Num dual var ", num_dualvar)
            println("Objective")
            dump(dual_objs)
            println("dual_A")
            println(dual_A)
            println("lhs,rhs idx")
            dump(gen.dual_ell_lhs_idxs)
            dump(gen.dual_ell_rhs_idxs)
            println("END DEBUG   :debug_printreform")
        end
    end  # end reformulation preparation
end

function generateCut(gen::GeneralOracle, master::Model, rm::Model, inds::Vector{Int}, active=false)
    # If not doing cuts...
    (!gen.use_cuts && !active) && return Any[]

    rd = getRobust(rm)
    master_sol = master.colVal
    new_cons = Any[]

    for con_ind in inds
        con = get_uncertain_constraint(rm, con_ind)

        # Update the cutting plane problem's objective, and solve
        cut_sense, unc_obj_coeffs, lhs_const = JuMPeR.build_cut_objective_sparse(con, master_sol)
        @setObjective(gen.cut_model, cut_sense, sum{u[2]*gen.cut_vars[u[1]], u=unc_obj_coeffs})
        cut_solve_status = solve(gen.cut_model, suppress_warnings=true)
        cut_solve_status != :Optimal &&
            error("GeneralOracle: cutting plane problem is infeasible or unbounded!")
        lhs_of_cut = getObjectiveValue(gen.cut_model) + lhs_const

        # SUBJECT TO CHANGE: active cut detection
        if active
            push!(rd.activecuts[con_ind], 
                cut_to_scen(gen.cut_model.colVal, 
                    check_cut_status(con, lhs_of_cut, gen.cut_tol) == :Active))
            continue
        end
        
        # Check violation
        if check_cut_status(con, lhs_of_cut, gen.cut_tol) != :Violate
            gen.debug_printcut && debug_printcut(rm,master,gen,lhs_of_cut,con,nothing)
            continue  # No violation, no new cut
        end
        
        # Create and add the new constraint
        new_con = JuMPeR.build_certain_constraint(master, con, gen.cut_model.colVal)
        gen.debug_printcut && debug_printcut(rm,master,gen,lhs_of_cut,con,new_con)
        push!(new_cons, new_con)
    end
    
    return new_cons
end


function debug_printcut(rm,m,w,lhs,con,new_con)
    println("BEGIN DEBUG :debug_printcut")
    #convert_model!(con, rm)
    println("  Constraint:  ", con)
    #convert_model!(con, m)
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




#############################################################################
# REFORMULATION
#############################################################################

function generateReform(gen::GeneralOracle, master::Model, rm::Model, inds::Vector{Int})
    # If not doing reform...
    gen.use_cuts && return 0
    # Apply the reformulation to all relevant constraints
    for ind in inds
        apply_reform(gen, master, rm, ind)
    end
    return length(inds)
end
    

function apply_reform(gen::GeneralOracle, master::Model, rm::Model, con_ind::Int)
    rd = getRobust(rm)
    con = get_uncertain_constraint(rm, con_ind)

    num_dualvar = gen.num_dualvar
    num_dualcon = rd.numUncs
    dual_A      = gen.dual_A
    dual_objs   = gen.dual_objs
    dual_vartype= gen.dual_vartype
    dual_contype= gen.dual_contype 
    # These are the objective coefficients of the cutting plane problem
    # that are the right-hand-side of the dual problem.
    dual_rhs    = [AffExpr() for i = 1:num_dualcon]
    # This is the reformulated "main" constraint that will replace the
    # existing uncertain constraint
    new_lhs     = AffExpr()

    # We will do all reformulation as if the constraint is a <=
    # constraint. This necessitates mulitplying through by -1
    # if the original constraint is not of this form
    sign_flip   = sense(con) == :(<=) ? +1.0 : -1.0
    # The right-hand-side of the new constraint
    new_rhs     = rhs(con) * sign_flip
    # The uncertain constraint left-hand-side
    orig_lhs    = con.terms

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
            rd.uncCat[unc] != :Cont && 
                error("No integer uncertainties allowed in reformulation!")
            push!(dual_rhs[unc], coeff * sign_flip, 
                                 Variable(master, var_col))
        end
    end
    # We also need the standalone (a^T u) not related to any variable
        for const_term_j = 1:length(orig_lhs.constant.coeffs)
            coeff = orig_lhs.constant.coeffs[const_term_j]
            unc   = orig_lhs.constant.vars[const_term_j].unc
            rd.uncCat[unc] != :Cont && 
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
        new_v = Variable(master,lbound,ubound,:Cont,var_name)
        push!(dual_vars, new_v)
        push!(new_lhs, dual_objs[ind], new_v)
    end

    gen.debug_printreform && println("DEBUG new_lhs ", new_lhs)

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
                      dual_vars[gen.dual_ell_rhs_idxs[ell_idx][el_matrix_row]] )
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
        beta_t = dual_vars[gen.dual_ell_lhs_idxs[ell_idx]]
        beta_zs = dual_vars[gen.dual_ell_rhs_idxs[ell_idx]]
        @addConstraint(master, sum{beta_z^2, beta_z in beta_zs} <= beta_t^2)
    end
    
    return true
end
