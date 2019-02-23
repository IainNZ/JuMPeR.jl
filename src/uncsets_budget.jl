#-----------------------------------------------------------------------
# JuMPeR  --  JuMP Extension for Robust Optimization
# http://github.com/IainNZ/JuMPeR.jl
#-----------------------------------------------------------------------
# Copyright (c) 2016: Iain Dunning
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
#-----------------------------------------------------------------------
# src/uncsets_budget.jl
# Defines the BudgetUncertaintySet, based on the set from the
# "Price of Robustness" paper by Bertsimas and Sim.
# Included by src/uncsets.jl
#-----------------------------------------------------------------------


"""
    BudgetUncertaintySet

Implements the "budget" uncertainty set from the 2004 paper "The Price
of Robustness" by Bertsimas and Sim. While this set is a polyhedron, it
is highly structured and admits an efficient cutting plane generation
method. The set is defined as follows: given uncertain parameters ξ,
nominal values μ, deviations σ, and a budget for number of allowed
deviations from nominal Γ, the budget uncertainty set can be written as

    U = { ξ | ξᵢ = μᵢ + σᵢzᵢ, ‖z‖₁ ≤ Γ, ‖z‖∞ ≤ 1 }

The data (Γ,μ,σ) is provided to the set when it is constructed. There
should be one value of μ and σ for every uncertain parameter in the problem.
Only the cutting plane method is implemented - there is nothing to gain
from a specific reformulation, as it is quite easy to achieve a similar effect
with the default `BasicUncertaintySet`:

    @uncertain(m,       ξ[1:n])
    @uncertain(m, -1 <= z[1:n] <= 1)  # ‖z‖∞ ≤ 1
    @constraint(m, norm(z, 1) <= Γ)  # ‖z‖₁ ≤ Γ
    for i in 1:n
        @constraint(m, ξ[i] == μ[i] + σ[i] * z[i])
    end

"""
mutable struct BudgetUncertaintySet <: AbstractUncertaintySet
    Γ::Int               # Budget of uncertainty
    μ::Vector{Float64}   # Nominal values for each parameter
    σ::Vector{Float64}   # Deviation values for each parameter
    ɛ::Float64           # Tolerance for new cutting planes
end
BudgetUncertaintySet(Γ, μ, σ) = BudgetUncertaintySet(Γ, μ, σ, 1e-4)


"""
    setup_set(BudgetUncertaintySet, ...)

BudgetUncertaintySet takes no action at this stage, but we must provide this
method as JuMPeR expects all uncertainty sets to implement it.
"""
setup_set(us::BudgetUncertaintySet, rm::Model, idxs::Vector{Int},
            scens_requested::Bool, other_prefs::Dict{Symbol,Any}) = nothing


"""
    get_worst_case_value(BudgetUncertaintySet, ...)

Internal function. For a given constraint and solution, find the worst-case
values of the uncertain parameters. This function is called by both
`generate_scenario` and `generate_cut`
"""
function get_worst_case_value(us::BudgetUncertaintySet, rm::Model, idx::Int)
    # Extract the RobustModelExt from the JuMP model - this contains all
    # the robust optimization-specific information
    rmext = get_robust(rm)::RobustModelExt
    # Get the UncConstraint object out of the model
    con = rmext.unc_constraints[idx]
    # Collect the xᵢ values for each uncertain parameter
    unc_x_vals = zeros(rmext.num_uncs)
    # We'll also calculate the "nominal" portion μᵢxᵢ
    nominal_value = 0.0
    # For every variable term in the constraint...
    for (unc_expr, var) in linearterms(con.terms)
        # Get the value of xᵢ in the current solution
        x_val = getvalue(var)
        # unc_expr is the coefficient on xᵢ. It might not be a single
        # uncertain parameter, e.g., (5a + 3b + 2)x, so we need to iterate
        # over this experession as well.
        for (coeff, unc) in linearterms(unc_expr)
            # Nominal portion first
            nominal_value += coeff * us.μ[unc.id] * x_val
            # Store the x value for this uncertain parameter
            unc_x_vals[unc.id] += coeff * x_val
        end
        # We may also have a deterministic part that factors only into the
        # nominal value calculation
        nominal_value += unc_expr.constant * x_val
    end
    # We may also have uncertain parameters not associated with any
    # variable. We can think of them as being associated with a variable
    # fixed at +1, so the code is similar to above. This often occurs due
    # to "right-hand-side" uncertain parameters.
    for (coeff, unc) in linearterms(con.terms.constant)
        nominal_value += coeff * us.μ[unc.id]
        unc_x_vals[unc.id] += coeff
    end
    nominal_value += con.terms.constant.constant
    # We now scale the x values by the deviations, and take the absolute
    # values - if σᵢ xᵢ is negative, we want to set ξᵢ=μᵢ-σᵢ, and if it
    # is positive, the opposite.
    scaled_vals = abs.(unc_x_vals) .* us.σ
    # We don't need to sort the list, just the permutation vector
    # of indices as if we had sorted. We then take the top Γ indices.
    max_inds = sortperm(scaled_vals)[(end - us.Γ + 1):end]
    # Determine the magnitude of the left-hand-side
    cut_value = nominal_value
    if JuMP.sense(con) == (:<=)
        # less-than-or-equal-to constraint - maximize left-hand-side
        cut_value += sum(scaled_vals[max_inds])
    else
        # greater-than-or-equal-to constraint - minimize left-hand-side
        cut_value -= sum(scaled_vals[max_inds])
    end
    # Calculate the uncertain parameter values
    uncvalues = copy(us.μ)
    for i in max_inds
        if unc_x_vals[i] > 0
            # Net positive coefficient on this uncertain parameter
            if JuMP.sense(con) == :(<=)
                uncvalues[i] += us.σ[i]  # Push up, LHS goes up
            else
                uncvalues[i] -= us.σ[i]  # Push down, LHS goes down
            end
        else
            # Net negative coefficient on this uncertain parameter
            if JuMP.sense(con) == :(<=)
                uncvalues[i] -= us.σ[i]  # Push down, LHS goes up
            else
                uncvalues[i] += us.σ[i]  # Push up, LHS goes down
            end
        end
    end
    # Return the LHS value (used by generate_cut) and the values of
    # the uncertain parameters (used by generate_cut and generate_scenario)
    return cut_value, uncvalues
end


"""
    generate_scenario(BudgetUncertaintySet, ...)

Wraps the results from `get_worst_case_value` in `Scenario` objects.
"""
function generate_scenario(us::BudgetUncertaintySet, rm::Model, idxs::Vector{Int})
    # We need to return one Scenario per constraint
    scens = Union{Scenario, Missing}[]
    for idx in idxs
        _, uncvalues = get_worst_case_value(us, rm, idx)
        push!(scens, Scenario(uncvalues))
    end
    return scens
end


"""
    generate_cut(BudgetUncertaintySet, ...)

Use the results from `get_worst_case_value` to determine if new new constraints
are needed, and return them if so.
"""
function generate_cut(us::BudgetUncertaintySet, rm::Model, idxs::Vector{Int})
    # Extract the RobustModelExt from the JuMP model - this contains all
    # the robust optimization-specific information
    rmext = get_robust(rm)::RobustModelExt
    # The vector of new constraints we will add
    new_cons = Any[]
    # For each constraint we need to generate a cut for
    for idx in idxs
        # Determine worst-case uncertain parameters
        cut_value, uncvalues = get_worst_case_value(us, rm, idx)
        # Get the UncConstraint object out of the robust model
        con = rmext.unc_constraints[idx]
        # Use a utility function from uncsets_util.jl to check the violation
        # that could be obtained by moving these uncertain parameters
        if check_cut_status(con, cut_value, us.ɛ) != :Violate
            # No violation, no new cut - cut maybe be :Active or :Slack
            continue  # try next constraint
        end
        # Build a deterministic constraint from the uncertain constraint
        new_con = build_certain_constraint(con, uncvalues)
        push!(new_cons, new_con)
    end
    return new_cons
end


"""
    generate_reform(BudgetUncertaintySet, ...)

Not implemented for this uncertainty set.
"""
generate_reform(us::BudgetUncertaintySet, rm::Model, idxs::Vector{Int}) = nothing
