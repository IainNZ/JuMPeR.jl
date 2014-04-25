#############################################################################
# JuMPeR
# Julia for Mathematical Programming - extension for Robust Optimization
# See http://github.com/IainNZ/JuMPeR.jl
#############################################################################
# BertSimOracle
# An oracle for the uncertainty set described in Bertsimas and Sim '04
# Although the set is polyhedral, it has special structure that is both
# tedious to manually construct but can also be exploited to efficiently
# find cuts.
#############################################################################

export BertSimOracle
type BertSimOracle <: AbstractOracle
    # The constraints associated with this oracle and the selected mode(s)
    # of operation for each constraint.
    cons::Vector{UncConstraint}
    con_modes::Vector{Dict{Symbol,Bool}}
    con_inds::Dict{Int,Int}
    any_reform::Bool
    any_cut::Bool
    setup_done::Bool


    Gamma::Int  # Only support integer values of Gamma for sake of cut
                # Could relax this, would make cutting plane more awkward
                # but is possible.
    means::Vector{Float64}
    devs::Vector{Float64}

    cut_tol
end
# Preferred constructor - just Gamma
BertSimOracle(Gamma::Int) =
    BertSimOracle(  UncConstraint[], Dict{Symbol,Bool}[], Dict{Int,Int}(),
                    false, false, false,
                    Gamma, Float64[], Float64[], 0.0)
# Default constructor - no uncertainty
BertSimOracle() = BertSimOracle(0)


# registerConstraint
# We must handle this constraint, and the users preferences have been
# communicated through prefs
function registerConstraint(w::BertSimOracle, con, ind::Int, prefs)
    w.cut_tol = get(prefs, :cut_tol, 1e-6)

    push!(w.cons, con)
    w.con_inds[ind] = length(w.cons)
    con_mode = Dict{Symbol,Bool}()
    #if :prefer_cuts in keys(prefs) && prefs[:prefer_cuts]
    #    con_mode = [:Cut => true, :Reform => false,
    #                :Sample => (:samples in keys(prefs)) ? 
    #                            prefs[:samples] > 0 : false  ]
        w.any_cut = true
    #else  # Default to reformulation
        #con_mode = [:Cut => false, :Reform => true, :Sample => false]
    #    w.any_reform = true
    #end
    con_mode = [:Cut => true, :Reform => true, :Sample => false]
    push!(w.con_modes, con_mode)
    return con_mode
end


# setup
# Analyze box on uncertainties to determine means/dev values for each uncertainty
function setup(w::BertSimOracle, rm::Model)
    rd = getRobust(rm)
    w.means = zeros(rd.numUncs)
    w.devs  = zeros(rd.numUncs)

    # Note, problematic we fail here - what if it doesn't appear in constraints
    # we are concerned with. Maybe this sohuld be done on a per-constraint basis?
    # Or at least not error out unless its tried to be used later on
    for i in 1:rd.numUncs
        if rd.uncUpper[i] == +Inf
            # No upper bound
            error("Uncertainty $(rm.uncNames[i]) (index $i) has no upper bound, cannot determine mean & deviation.")
        elseif rd.uncLower[i] == -Inf
            error("Uncertainty $(rm.uncNames[i]) (index $i) has no lower bound, cannot determine mean & deviation.")
        end
        w.means[i] = (rd.uncUpper[i] + rd.uncLower[i])/2
        w.devs[i]  = (rd.uncUpper[i] - rd.uncLower[i])/2
    end
end


function generateCut(w::BertSimOracle, rm::Model, ind::Int, m::Model, cb=nothing)

    # If not doing cuts for this one, just skip
    con_ind = w.con_inds[ind]
    if !w.con_modes[con_ind][:Cut]
        return 0
    end

    # Given a constraint
    # u_1 x_1 + .... u_n x_n <= b
    # with this uncertainty set, it is sufficient to:
    # 1. Calculate abs(x_i) * dev_i for each i
    # 2. Sort ascending
    # 3. If constraint can be violated by setting Gamma at bounds
    #    add cut.
    # TODO: Relax assumption of one uncertain per variable?
    master_sol = m.colVal

    absx_devs = Float64[]       # kth term is x[j]*devs[i]
    uncx_inds = Int[]           # k

    con = w.cons[con_ind]
    orig_lhs = con.terms
    nom_val  = 0.0
    for var_ind = 1:length(orig_lhs.vars)
        col      = orig_lhs.vars[var_ind].col
        num_uncs = length(orig_lhs.coeffs[var_ind].vars)
        coeff_val = 0.0
        if num_uncs > 1
            # More than one uncertain on this variable
            # Not supported
            error("BertSimOracle only supports one uncertain coefficient per variable")
        elseif num_uncs == 1
            unc = orig_lhs.coeffs[var_ind].vars[1].unc
            push!(absx_devs, abs(master_sol[col]) * w.devs[unc])
            push!(uncx_inds, var_ind)
            coeff_val += orig_lhs.coeffs[var_ind].coeffs[1] * w.means[unc]
        end
        coeff_val += orig_lhs.coeffs[var_ind].constant
        nom_val += coeff_val * master_sol[col]
    end
    if length(orig_lhs.constant.vars) >= 1
        error("BertSimOracle doesn't support uncertain not attached to variable yet")
    end
    nom_val += orig_lhs.constant.constant

    # Now we have the lists, we can sort them in order
    perm = sortperm(absx_devs)  # Ascending
    max_inds = uncx_inds[perm[(end-w.Gamma+1):end]]
    
    # Check violation
    cut_val = nom_val + 
                ((sense(con) == :<=) ? +1.0 : -1.0) * sum(absx_devs[max_inds])
    
    if !is_constraint_violated(con, cut_val, w.cut_tol)
        return 0  # No violation, no new cut
    end

    # Add new constraint
    unc_val = w.means[:]
    for p in 1:length(perm)
        var_ind = uncx_inds[p]
        unc     = orig_lhs.coeffs[var_ind].vars[1].unc
        col     = orig_lhs.vars[var_ind].col
        # Whether we add or remove a deviation depends on both the 
        # constraint sense and the sign of x
        sign    = sense(con) == :(<=) ? (master_sol[col] >= 0 ? +1.0 : -1.0) :
                                        (master_sol[col] >= 0 ? -1.0 : +1.0)
        unc_val[unc] += sign * w.devs[unc]
    end
    new_con = JuMPeR.build_certain_constraint(con, unc_val)
    cb == nothing ? addConstraint(m, new_con) :
                    addLazyConstraint(cb, new_con)

    return 1
end


function generateReform(w::BertSimOracle, rm::Model, ind::Int, m::Model)
    # If not doing reform for this one, just skip
    con_ind = w.con_inds[ind]
    if !w.con_modes[con_ind][:Reform]
        return false
    end
    
    # Generate an initial cut for each constraint using the nominal values
    # To help make initial cut round more useful
    con = w.cons[con_ind]
        orig_lhs = con.terms
        new_lhs = AffExpr(orig_lhs.vars,
                          [orig_lhs.coeffs[i].constant for i in 1:length(orig_lhs.vars)],
                          orig_lhs.constant.constant)
        for var_ind = 1:length(orig_lhs.vars)
            col      = orig_lhs.vars[var_ind].col
            num_uncs = length(orig_lhs.coeffs[var_ind].vars)
            if num_uncs > 1
                # More than one uncertain on this variable - not supported
                error("BertSimOracle only supports one uncertain coefficient per variable")
            elseif num_uncs == 1
                unc   = orig_lhs.coeffs[var_ind].vars[1].unc
                coeff = orig_lhs.coeffs[var_ind].coeffs[1]
                new_lhs.coeffs[var_ind] += coeff * w.means[unc]
            end
        end
        if sense(con) == :<=
            @addConstraint(m, new_lhs <= con.ub)
        else
            @addConstraint(m, new_lhs >= con.lb)
        end
    return true
end