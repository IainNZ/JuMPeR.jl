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

    # Mappings from full set to components
    unc_to_comp::Vector{Int}
    unc_to_comp_unc::Vector{Vector{Int}}
    con_to_comp::Vector{Int}

    # Options
    debug_printcut::Bool
end
# Default constructor
GeneralGraphOracle() = 
    GeneralGraphOracle( Model[], Vector{Variable}[], 0.0,
                        Int[], Vector{Int}[], Int[],
                        false)


# registerConstraint
# We must handle this constraint, and the users preferences have been
# communicated through prefs. We don't need to take any action here.
registerConstraint(gen::GeneralGraphOracle, rm::Model, ind::Int, prefs) = nothing


# setup
# Generate the cutting plane model or precompute the reformulation structure.
function setup(gen::GeneralGraphOracle, rm::Model, prefs)

    # Extract preferences we care about
    gen.cut_tol           = get(prefs, :cut_tol, 1e-6)
    gen.debug_printcut    = get(prefs, :debug_printcut, false)

    rd = getRobust(rm)

    # Analyze components
    gen.unc_to_comp, gen.con_to_comp = detect_components(rd.numUncs, rd.uncertaintyset)
    num_components = maximum(gen.unc_to_comp)
    gen.debug_printcut && println("GeneralGraphOracle: $num_components components detected.")
    comp_to_unc = [Int[] for c in 1:num_components]
    for i = 1:rd.numUncs
        push!(comp_to_unc[gen.unc_to_comp[i]], i)
    end

    # Cutting plane setup
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
        m.colCat   = rd.uncCat  [comp_to_unc[comp]]
        cut_vars   = [Variable(m, i) for i = 1:length(comp_to_unc[comp])]
        # Polyhedral constraints only right now
        for con_ind in 1:length(rd.uncertaintyset)
            gen.con_to_comp[con_ind] != comp && continue
            c = rd.uncertaintyset[con_ind]
            newcon = LinearConstraint(AffExpr(), c.lb, c.ub)
            newcon.terms.coeffs = c.terms.coeffs
            newcon.terms.vars   = [cut_vars[unc_to_comp_unc[u.unc]] for u in c.terms.vars]
            push!(m.linconstr, newcon)
        end

        push!(gen.unc_to_comp_unc, unc_to_comp_unc)
        push!(gen.comp_cut_model, m)
        push!(gen.comp_cut_vars, cut_vars)
    end
end


function generateReform(gen::GeneralGraphOracle, master::Model, rm::Model, inds::Vector{Int})
    return 0
end


function generateCut(gen::GeneralGraphOracle, master::Model, rm::Model, inds::Vector{Int}, active=false)

    rd = getRobust(rm)
    master_sol = master.colVal
    new_cons = {}

    for con_ind in inds
        con = get_uncertain_constraint(rm, con_ind)

        # Pull out right model for this cut
        comp = 0
        unc_lhs = con.terms
        for var_ind = 1:length(unc_lhs.vars)
            uaff = unc_lhs.coeffs[var_ind]
            for unc_ind = 1:length(uaff.vars)
                comp = gen.unc_to_comp[uaff.vars[unc_ind].unc]
                break
            end
        end
        # Uncertains not attached to variables
        if comp == 0
            uaff = unc_lhs.constant
            for unc_ind = 1:length(uaff.vars)
                comp = gen.unc_to_comp[uaff.vars[unc_ind].unc]
                break
            end
        end

        @assert comp != 0
        cut_model       = gen.comp_cut_model[comp]
        cut_vars        = gen.comp_cut_vars [comp]
        unc_to_comp_unc = gen.unc_to_comp_unc[comp]

        # Update the cutting plane problem's objective, and solve
        cut_sense, unc_obj_coeffs, lhs_const = JuMPeR.build_cut_objective_sparse(con, master_sol)
        
        cut_obj = AffExpr()
        for uoc in unc_obj_coeffs
            gen.unc_to_comp[uoc[1]] != comp && continue
            push!(cut_obj, uoc[2], cut_vars[unc_to_comp_unc[uoc[1]]])
        end
        @setObjective(cut_model, cut_sense, cut_obj)

        cut_solve_status = solve(cut_model, suppress_warnings=true)
        cut_solve_status != :Optimal && error("Cutting plane problem infeasible or unbounded!")
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
            #gen.debug_printcut && debug_printcut(rm,master,gen,lhs_of_cut,con,nothing)
            continue  # No violation, no new cut
        end
        
        # Create and add the new constraint
        new_con = JuMPeR.build_certain_constraint(master, con, full_colVal)
        #gen.debug_printcut && debug_printcut(rm,master,gen,lhs_of_cut,con,new_con)
        push!(new_cons, new_con)
    end
    
    return new_cons
end