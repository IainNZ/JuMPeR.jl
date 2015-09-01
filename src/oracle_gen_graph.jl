#############################################################################
# JuMPeR
# Julia for Mathematical Programming - extension for Robust Optimization
# See http://github.com/IainNZ/JuMPeR.jl
#############################################################################
# GeneralGraphOracle
# The GeneralOracle but with extra smarts to split up the uncertainty set
# if it is actually multiple uncertainty sets. Doesn't yet support
# reformulations but could in theory - just a matter of coding it.
#############################################################################

type GeneralGraphOracle <: AbstractOracle
    # Cutting plane models
    comp_cut_model::Vector{Model}
    comp_cut_vars::Vector{Vector{Variable}}
    cut_tol::Float64
    use_cuts::Bool

    # Reformulation (see setup() for details)
    comp_num_dualvar::Vector{Int}
    comp_num_dualcon::Vector{Int}
    comp_dual_A::Vector{Vector{Vector{@compat Tuple{Int,Float64}}}}
    comp_dual_objs::Vector{Vector{Float64}}
    comp_dual_vartype::Vector{Vector{Symbol}}
    comp_dual_contype::Vector{Vector{Symbol}}

    # Mappings from full set to components
    unc_to_comp::Vector{Int}
    unc_to_comp_unc::Vector{Vector{Int}}
    con_to_comp::Vector{Int}
end
# Default constructor
GeneralGraphOracle() = 
    GeneralGraphOracle( Model[], Vector{Variable}[], 0.0, false,
                        Int[], Int[], Vector{Vector{@compat Tuple{Int,Float64}}}[], Vector{Float64}[],
                        Vector{Symbol}[], Vector{Symbol}[],
                        Int[], Vector{Int}[], Int[])


# registerConstraint
# We must handle this constraint, and the users preferences have been
# communicated through prefs. We don't need to take any action here.
registerConstraint(gen::GeneralGraphOracle, rm::Model, ind::Int, prefs) = nothing


# setup
# Generate the cutting plane model or precompute the reformulation structure.
function setup(gen::GeneralGraphOracle, rm::Model, prefs)

    # Extract preferences we care about
    gen.use_cuts          = get(prefs, :prefer_cuts, false)
    gen.cut_tol           = get(prefs, :cut_tol, 1e-6)

    rd = getRobust(rm)

    # Analyze components
    gen.unc_to_comp, gen.con_to_comp = detect_components(rd.numUncs, rd.uncertaintyset)
    num_components = maximum(gen.unc_to_comp)
    #@show num_components
    #gen.debug_printcut && println("GeneralGraphOracle: $num_components components detected.")
    comp_to_unc = [Int[] for c in 1:num_components]
    for i = 1:rd.numUncs
        push!(comp_to_unc[gen.unc_to_comp[i]], i)
    end

    # Cutting plane setup
    if gen.use_cuts || prefs[:active_cuts]
        for comp = 1:num_components
            m = Model(solver=rd.cutsolver == nothing ? rm.solver : rd.cutsolver)
            unc_to_comp_unc = zeros(Int, rd.numUncs)
            pos = 0
            for i = 1:rd.numUncs
                gen.unc_to_comp[i] != comp && continue
                pos += 1
                unc_to_comp_unc[i] = pos
            end
            m.numCols  = length(comp_to_unc[comp])
            m.colNames = rd.uncNames[comp_to_unc[comp]]
            m.colLower = rd.uncLower[comp_to_unc[comp]]
            m.colUpper = rd.uncUpper[comp_to_unc[comp]]
            m.colCat   = rd.uncCat[  comp_to_unc[comp]]
            cut_vars   = [Variable(m, i) for i = 1:length(comp_to_unc[comp])]
            # Polyhedral constraints only right now
            for con_ind in 1:length(rd.uncertaintyset)
                gen.con_to_comp[con_ind] != comp && continue
                c = rd.uncertaintyset[con_ind]
                newcon = LinearConstraint(AffExpr(), c.lb, c.ub)
                newcon.terms.coeffs = c.terms.coeffs
                newcon.terms.vars   = [cut_vars[unc_to_comp_unc[u.id]] for u in c.terms.vars]
                push!(m.linconstr, newcon)
            end

            push!(gen.unc_to_comp_unc, unc_to_comp_unc)
            push!(gen.comp_cut_model, m)
            push!(gen.comp_cut_vars, cut_vars)
        end
    end

    # Reformulation setup
    if !gen.use_cuts
        for comp = 1:num_components
            # Get constraints for this component
            comp_uncset = Any[]
            for ci = 1:length(rd.uncertaintyset)
                if gen.con_to_comp[ci] == comp
                    push!(comp_uncset, rd.uncertaintyset[ci])
                end
            end
            # Get uncertains for this component
            unc_to_comp_unc = zeros(Int, rd.numUncs)
            pos = 0
            for i = 1:rd.numUncs
                gen.unc_to_comp[i] != comp && continue
                pos += 1
                unc_to_comp_unc[i] = pos
            end

            num_dualvar = length(comp_uncset)       # WAS: length(rd.uncertaintyset)
            num_dualcon = pos                       # WAS: rd.numUncs
            dual_A      = [(@compat Tuple{Int,Float64})[] for i = 1:num_dualcon]   # WAS: rd.numUncs
            dual_objs   = zeros(num_dualvar)
            dual_vartype = [:>= for i = 1:num_dualvar]
            dual_contype = [:>= for i = 1:num_dualcon]

            for uncset_i = 1:length(comp_uncset)    # WAS: length(rd.uncertaintyset)
                lhs = comp_uncset[uncset_i].terms   # WAS: rd.uncertaintyset[uncset_i].terms
                for unc_j = 1:length(lhs.vars)
                    push!(dual_A[unc_to_comp_unc[lhs.vars[unc_j].id]],  # WAS: [lhs.vars[unc_j].unc],
                                (uncset_i, lhs.coeffs[unc_j]))
                end
                dual_objs[uncset_i]    = rhs(comp_uncset[uncset_i])
                dual_vartype[uncset_i] = sense(comp_uncset[uncset_i])
            end

            for unc_i = 1:rd.numUncs
                gen.unc_to_comp[unc_i] != comp && continue
                # If it has a lower bound...
                if rd.uncLower[unc_i] != -Inf
                    num_dualvar += 1
                    push!(dual_A[unc_to_comp_unc[unc_i]], (num_dualvar, 1.0))
                    push!(dual_objs,     rd.uncLower[unc_i])
                    push!(dual_vartype,  :(>=))
                end
                # If it has an upper bound
                if rd.uncUpper[unc_i] != +Inf
                    num_dualvar += 1
                    push!(dual_A[unc_to_comp_unc[unc_i]], (num_dualvar, 1.0))
                    push!(dual_objs,     rd.uncUpper[unc_i])
                    push!(dual_vartype,  :(<=))
                end
                # Now we just treat the variable as being free
                dual_contype[unc_to_comp_unc[unc_i]] = :(==)
            end

            #=println("BEGIN DEBUG :debug_printreform")
            println("Component: ", comp)
            println("Mapping: ", unc_to_comp_unc)
            println("Num dual var ", num_dualvar)
            println("Objective")
            dump(dual_objs)
            println("dual_A")
            println(dual_A)
            println("END DEBUG   :debug_printreform")=#

            push!(gen.comp_num_dualvar  , num_dualvar   )
            push!(gen.comp_num_dualcon  , num_dualcon   )
            push!(gen.comp_dual_A       , dual_A        )
            push!(gen.comp_dual_objs    , dual_objs     )
            push!(gen.comp_dual_vartype , dual_vartype  )
            push!(gen.comp_dual_contype , dual_contype  )
        end # next component
    end
end


function generateReform(gen::GeneralGraphOracle, master::Model, rm::Model, inds::Vector{Int})
    gen.use_cuts && return 0
    rd = getRobust(rm)
    for con_ind in inds
        con = get_uncertain_constraint(rm, con_ind)

        # Pull out right component
        comp = 0
        unc_lhs = con.terms
        for var_ind = 1:length(unc_lhs.vars)
            uaff = unc_lhs.coeffs[var_ind]
            for unc_ind = 1:length(uaff.vars)
                comp = gen.unc_to_comp[uaff.vars[unc_ind].id]
                break
            end
        end
        if comp == 0
            uaff = unc_lhs.constant
            for unc_ind = 1:length(uaff.vars)
                comp = gen.unc_to_comp[uaff.vars[unc_ind].id]
                break
            end
        end
        #println("genReform $comp")

        # Create mapping again
        unc_to_comp_unc = zeros(Int, rd.numUncs)
        pos = 0
        for i = 1:rd.numUncs
            gen.unc_to_comp[i] != comp && continue
            pos += 1
            unc_to_comp_unc[i] = pos
        end

        num_dualvar = gen.comp_num_dualvar[comp]
        num_dualcon = gen.comp_num_dualcon[comp]
        dual_A      = gen.comp_dual_A[comp]
        dual_objs   = gen.comp_dual_objs[comp]
        dual_vartype= gen.comp_dual_vartype[comp]
        dual_contype= gen.comp_dual_contype[comp]

        dual_rhs    = [AffExpr() for i = 1:num_dualcon]
        new_lhs     = AffExpr()
        sign_flip   = sense(con) == :(<=) ? +1.0 : -1.0
        new_rhs     = rhs(con) * sign_flip
        orig_lhs    = con.terms

        # Certain terms
        num_lhs_terms = length(orig_lhs.vars)
        for term_i = 1:num_lhs_terms
            var_col = orig_lhs.vars[term_i].col
            if orig_lhs.coeffs[term_i].constant != 0.0
                push!(new_lhs,  orig_lhs.coeffs[term_i].constant * sign_flip,
                                Variable(master, var_col) )
            end
        end

        # RHS from var coefficients
        for term_i = 1:num_lhs_terms
            var_col     = orig_lhs.vars[term_i].col
            term_coeff  = orig_lhs.coeffs[term_i]
            for coeff_term_j = 1:length(term_coeff.coeffs)
                coeff = term_coeff.coeffs[coeff_term_j]
                unc   = term_coeff.vars[coeff_term_j].id
                push!(dual_rhs[unc_to_comp_unc[unc]], 
                                    coeff * sign_flip, 
                                     Variable(master, var_col))
            end
        end
        # From constant term
            for const_term_j = 1:length(orig_lhs.constant.coeffs)
                coeff = orig_lhs.constant.coeffs[const_term_j]
                unc   = orig_lhs.constant.vars[const_term_j].id
                dual_rhs[unc_to_comp_unc[unc]].constant += coeff * sign_flip
            end

        # Make dual vars
        dual_vars = Variable[]
        for ind = 1:num_dualvar
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


        @addConstraint(master, new_lhs <= new_rhs)

        for unc = 1:num_dualcon
            new_lhs = AffExpr()
            # Add the linear constraint segment
            for pair in dual_A[unc]
                push!(new_lhs, pair[2], dual_vars[pair[1]])
            end

            dualtype = dual_contype[unc]
            if     dualtype == :(==)
                @addConstraint(master, new_lhs == dual_rhs[unc])
            elseif dualtype == :(<=)
                @addConstraint(master, new_lhs <= dual_rhs[unc])
            else
                @addConstraint(master, new_lhs >= dual_rhs[unc])
            end
        end
    end
    return length(inds)
end


function generateCut(gen::GeneralGraphOracle, master::Model, rm::Model, inds::Vector{Int}, active=false)
    # If not doing cuts...
    (!gen.use_cuts && !active) && return Any[]

    rd = getRobust(rm)
    master_sol = master.colVal
    new_cons = Any[]

    for con_ind in inds
        con = get_uncertain_constraint(rm, con_ind)

        # Pull out right model for this cut
        comp = 0
        unc_lhs = con.terms
        for var_ind = 1:length(unc_lhs.vars)
            uaff = unc_lhs.coeffs[var_ind]
            for unc_ind = 1:length(uaff.vars)
                comp = gen.unc_to_comp[uaff.vars[unc_ind].id]
                break
            end
        end
        # Uncertains not attached to variables
        if comp == 0
            uaff = unc_lhs.constant
            for unc_ind = 1:length(uaff.vars)
                comp = gen.unc_to_comp[uaff.vars[unc_ind].id]
                break
            end
        end

        @assert comp != 0
        cut_model       = gen.comp_cut_model[ comp]
        cut_vars        = gen.comp_cut_vars[  comp]
        unc_to_comp_unc = gen.unc_to_comp_unc[comp]

        # Update the cutting plane problem's objective, and solve
        cut_sense, unc_obj_coeffs, lhs_const = JuMPeR.build_cut_objective_sparse(con, master_sol)
        
        cut_obj = AffExpr()
        for uoc in unc_obj_coeffs
            gen.unc_to_comp[uoc[1]] != comp && continue
            push!(cut_obj, uoc[2], cut_vars[unc_to_comp_unc[uoc[1]]])
        end
        @setObjective(cut_model, cut_sense, cut_obj)

        cut_solve_status = solve(cut_model) #, suppress_warnings=true)
        cut_solve_status != :Optimal &&
            error("GeneralGraphOracle: Cutting plane problem infeasible or unbounded!")
        lhs_of_cut = getObjectiveValue(cut_model) + lhs_const
        
        # Expand cut solution to full uncertainty set space
        full_colVal = zeros(rd.numUncs)
        for i in 1:rd.numUncs
            gen.unc_to_comp[i] != comp && continue
            full_colVal[i] = cut_model.colVal[unc_to_comp_unc[i]]
        end

        # SUBJECT TO CHANGE: active cut detection
        if active
            push!(rd.activecuts[con_ind], 
                cut_to_scen(full_colVal, 
                    check_cut_status(con, lhs_of_cut, gen.cut_tol) == :Active))
            continue
        end
        
        # Check violation
        if check_cut_status(con, lhs_of_cut, gen.cut_tol) != :Violate
            continue  # No violation, no new cut
        end
        
        # Create and add the new constraint
        new_con = JuMPeR.build_certain_constraint(master, con, full_colVal)
        push!(new_cons, new_con)
    end
    
    return new_cons
end