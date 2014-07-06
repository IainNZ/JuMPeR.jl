#############################################################################
# JuMPeR
# Julia for Mathematical Programming - extension for Robust Optimization
# See http://github.com/IainNZ/JuMPeR.jl
#############################################################################
# GeneralGraphOracle
# EXPERIMENTAL
# The default oracle - uses uncertainty bounds and the uncertainty set
# built up with linear constraints and ellipses and either use cutting planes
# or reformulate to obtain a deterministic problem.
#############################################################################

type GeneralGraphOracle <: AbstractOracle
    # The constraints associated with this oracle and the selected mode(s)
    # of operation for each constraint.
    cons::Vector{UncConstraint}
    con_modes::Vector{Dict{Symbol,Bool}}
    con_inds::Dict{Int,Int}
    setup_done::Bool

    # Multiple cutting plane models
    unc_to_comp::Vector{Int}
    unc_to_comp_unc::Vector{Vector{Int}}
    con_to_comp::Vector{Int}
    comp_cut_model::Vector{Model}
    comp_cut_vars::Vector{Vector{Variable}}

    # Other options
    debug_printcut::Bool
    cut_tol::Float64
end
# Default constructor
GeneralGraphOracle() = GeneralGraphOracle(
                        UncConstraint[], Dict{Symbol,Bool}[], Dict{Int,Int}(), false,
                        Int[], Vector{Int}[], Int[], Model[], Vector{Variable}[],
                        false, 1e-6)


# registerConstraint
# We must handle this constraint, and the users preferences have been
# communicated through prefs
function registerConstraint(w::GeneralGraphOracle, con, ind::Int, prefs)
    con_mode = [:Cut    =>  get(prefs, :prefer_cuts, false), 
                :Reform => !get(prefs, :prefer_cuts, false)]
    push!(w.con_modes, con_mode)
    push!(w.cons, con)
    w.con_inds[ind] = length(w.cons)

    # Extract preferences we care about
    w.debug_printcut    = get(prefs, :debug_printcut, false)
    w.cut_tol           = get(prefs, :cut_tol, 1e-6)

    return con_mode
end


# setup
# Now that all modes of operation have been selected, generate the cutting
# plane model and/or the reformulation
function setup(w::GeneralGraphOracle, rm::Model)
    w.setup_done && return
    rd = getRobust(rm)
    any_cut = any_reform = false
    for con_mode in w.con_modes
        any_cut    |= con_mode[:Cut]
        any_reform |= con_mode[:Reform]
    end

    # Analyze components
    w.unc_to_comp, w.con_to_comp = detect_components(rd.numUncs, rd.uncertaintyset)
    num_components = maximum(w.unc_to_comp)
    println("GeneralGraphOracle: $num_components components detected.")
    comp_to_unc = [Int[] for c in 1:num_components]
    for i = 1:rd.numUncs
        push!(comp_to_unc[w.unc_to_comp[i]], i)
    end

    for comp = 1:num_components
        m = Model(solver=rd.cutsolver == nothing ? rm.solver : rd.cutsolver)
        unc_to_comp_unc = zeros(Int, rd.numUncs)
        pos = 0
        for i = 1:rd.numUncs
            w.unc_to_comp[i] != comp && continue
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
            w.con_to_comp[con_ind] != comp && continue
            c = rd.uncertaintyset[con_ind]
            newcon = LinearConstraint(AffExpr(), c.lb, c.ub)
            newcon.terms.coeffs = c.terms.coeffs
            newcon.terms.vars   = [cut_vars[unc_to_comp_unc[u.unc]] for u in c.terms.vars]
            push!(m.linconstr, newcon)
        end
        #println("COMP $comp")
        #println(m)
        #println(unc_to_comp_unc)

        push!(w.unc_to_comp_unc, unc_to_comp_unc)
        push!(w.comp_cut_model, m)
        push!(w.comp_cut_vars, cut_vars)
    end

    w.setup_done = true
end


function generateCut(w::GeneralGraphOracle, rm::Model, ind::Int, m::Model, cb=nothing, active=false)
    # If not doing cuts for this one, just skip
    con_ind = w.con_inds[ind]
    if !w.con_modes[con_ind][:Cut] && !active
        return 0
    end
    con = w.cons[con_ind]
    rd = getRobust(rm)
    
    master_sol = m.colVal

    # Pull out right model for this cut
    comp = 0
    unc_lhs = con.terms
    for var_ind = 1:length(unc_lhs.vars)
        uaff = unc_lhs.coeffs[var_ind]
        for unc_ind = 1:length(uaff.vars)
            comp = w.unc_to_comp[uaff.vars[unc_ind].unc]
            break
        end
    end
    # Uncertains not attached to variables
    if comp == 0
        uaff = unc_lhs.constant
        for unc_ind = 1:length(uaff.vars)
            comp = w.unc_to_comp[uaff.vars[unc_ind].unc]
            break
        end
    end

    @assert comp != 0
    cut_model       = w.comp_cut_model[comp]
    cut_vars        = w.comp_cut_vars [comp]
    unc_to_comp_unc = w.unc_to_comp_unc[comp]

    # Update the cutting plane problem's objective, and solve
    cut_sense, unc_obj_coeffs, lhs_const = JuMPeR.build_cut_objective(con, master_sol)

    cut_obj = AffExpr()
    for uoc in unc_obj_coeffs
        w.unc_to_comp[uoc[1]] != comp && continue
        push!(cut_obj, uoc[2], cut_vars[unc_to_comp_unc[uoc[1]]])
    end
    @setObjective(cut_model, cut_sense, cut_obj)

    cut_solve_status = solve(cut_model, suppress_warnings=true)
    cut_solve_status != :Optimal && error("Cutting plane problem infeasible or unbounded!")
    lhs_of_cut = getObjectiveValue(cut_model) + lhs_const
    full_colVal = zeros(rd.numUncs)
    for i in 1:rd.numUncs
        w.unc_to_comp[i] != comp && continue
        full_colVal[i] = cut_model.colVal[unc_to_comp_unc[i]]
    end

    # TEMPORARY: active cut detection
    if active        
        push!(rd.activecuts[ind], 
            cut_to_scen(full_colVal, 
                check_cut_status(con, lhs_of_cut, w.cut_tol) == :Active))
    end
    
    # Check violation
    if check_cut_status(con, lhs_of_cut, w.cut_tol) != :Violate
        return 0  # No violation, no new cut
    end
    
    # Create and add the new constraint
    new_con = JuMPeR.build_certain_constraint(con, full_colVal)
    cb == nothing ? addConstraint(m, new_con) :
                    addLazyConstraint(cb, new_con)
    return 1
end


function generateReform(w::GeneralGraphOracle, rm::Model, ind::Int, master::Model)
    return false
end
