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
    Gamma::Int              # Only support integer values of Gamma (simpler)
    noms::Vector{Float64}   # Nominal values for each uncertainty
    devs::Vector{Float64}   # Deviation values for each uncertainty
    cut_tol::Float64
end
# Preferred constructor - just Gamma
BertSimOracle(Gamma::Int) = BertSimOracle(Gamma, Float64[], Float64[], 1e-6)
# Default constructor - no uncertainty
BertSimOracle() = BertSimOracle(0)


# registerConstraint
# We need to validate that every uncertain that appears in this 
# constraint has bound information. We'll also calculate the nominal
# and deviation while we are at it
function registerConstraint(bs::BertSimOracle, rm::Model, ind::Int, prefs)
    rd = getRobust(rm)

    # If we haven't allocated space for the nominals/deviations yet,
    # do so now
    if length(bs.noms) == 0
        bs.noms = fill(NaN, rd.numUncs)
        bs.devs = fill(NaN, rd.numUncs)
    end

    # Check every uncertain in this constraint to ensure it has bounds
    # and if so, calculate the mean and deviation values
    con = get_uncertain_constraint(rm, ind)
    for term_ind in 1:length(con.terms.vars)
        for unc_ind in 1:length(con.terms.coeffs[term_ind].vars)
            unc = con.terms.coeffs[term_ind].vars[unc_ind].unc
            if isnan(bs.noms[unc])
                lower, upper = rd.uncLower[unc], rd.uncUpper[unc]
                upper == +Inf &&
                    error("Unc. $(rm.uncNames[i]) has no upper bound, cannot determine nominal & deviation.")
                lower == -Inf &&
                    error("Unc. $(rm.uncNames[i]) has no lower bound, cannot determine nominal & deviation.")
                bs.noms[unc] = (upper + lower)/2
                bs.devs[unc] = (upper - lower)/2
            end
        end
    end
end


# setup
# Most work has been done in registerConstraint
function setup(bs::BertSimOracle, rm::Model, prefs)
    bs.cut_tol = get(prefs, :cut_tol, 1e-6)
end


# generateCut
# Generate cuts by doing a sort based on current value of the master problem
# Given a constraint
# u_1 x_1 + .... u_n x_n <= b
# with this uncertainty set, it is sufficient to:
# 1. Calculate abs(x_i) * dev_i for each i
# 2. Sort ascending
# 3. If constraint can be violated by setting Gamma at bounds, add cut.
# TODO: Relax assumption of one uncertain per variable?
function generateCut(bs::BertSimOracle, master::Model, rm::Model, inds::Vector{Int}, active=false)
    master_sol = master.colVal
    new_cons = {}

    for con_ind in inds
        absx_devs = Float64[]       # kth term is x[j]*devs[i]
        uncx_inds = Int[]           # k

        con = get_uncertain_constraint(rm, con_ind)
        orig_lhs = con.terms
        nom_val  = 0.0

        # Calculate the nominal value of the LHS constraint, and how much we
        # can move that LHS value by adjusting the uncertains.
        # For every variable in the constraint...
        for var_ind = 1:length(orig_lhs.vars)
            num_uncs = length(orig_lhs.coeffs[var_ind].vars)
            # Only support one uncertain per variable
            num_uncs > 1 && error("BertSimOracle only supports one uncertain per variable")
            
            col = orig_lhs.vars[var_ind].col
            nom_coeff_val = orig_lhs.coeffs[var_ind].constant
            if num_uncs == 1
                unc = orig_lhs.coeffs[var_ind].vars[1].unc
                # Store |x|*deviation
                push!(absx_devs, abs(master_sol[col]) * bs.devs[unc])
                push!(uncx_inds, var_ind)
                nom_coeff_val += orig_lhs.coeffs[var_ind].coeffs[1] * bs.noms[unc]
            end
            nom_val += nom_coeff_val * master_sol[col]
        end
        length(orig_lhs.constant.vars) >= 1 &&
            error("BertSimOracle doesn't support unattached uncertainties.")
        nom_val += orig_lhs.constant.constant

        # Sort the list of how much we can change LHS in ascending order
        perm = sortperm(absx_devs)
        # Obtain the top Gamma indices 
        max_inds = uncx_inds[perm[(end-bs.Gamma+1):end]]
        
        # Check violation that would be obtained from moving these in the
        # adversarial direction
        cut_val = nom_val + 
                    ((sense(con) == :<=) ? +1.0 : -1.0) * sum(absx_devs[max_inds])
        if check_cut_status(con, cut_val, bs.cut_tol) != :Violate
            # No violation, no new cut
            continue
        end

        # Violation -> add new constraint
        unc_val = bs.noms[:]
        for p in 1:length(perm)
            var_ind = uncx_inds[p]
            unc     = orig_lhs.coeffs[var_ind].vars[1].unc
            col     = orig_lhs.vars[var_ind].col
            # Whether we add or remove a deviation depends on both the 
            # constraint sense and the sign of x
            sign    = sense(con) == :(<=) ? (master_sol[col] >= 0 ? +1.0 : -1.0) :
                                            (master_sol[col] >= 0 ? -1.0 : +1.0)
            unc_val[unc] += sign * bs.devs[unc]
        end
        new_con = JuMPeR.build_certain_constraint(con, unc_val)
        push!(new_cons, new_con)
    end

    return new_cons
end


# generateReform
# Not implemented yet for this oracle
function generateReform(bs::BertSimOracle, master::Model, rm::Model, inds::Vector{Int})
    return 0 
end