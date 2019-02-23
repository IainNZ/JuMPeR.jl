#-----------------------------------------------------------------------
# JuMPeR  --  JuMP Extension for Robust Optimization
# http://github.com/IainNZ/JuMPeR.jl
#-----------------------------------------------------------------------
# Copyright (c) 2016: Iain Dunning
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
#-----------------------------------------------------------------------
# src/uncsets.jl
# Defines the AbstractUncertaintySet interface, which controls how
# constraints with uncertain parameters are addressed when solving a
# robust optimization problem.
# Included by src/JuMPeR.jl
#-----------------------------------------------------------------------

# AbstractUncertaintySet is defined in src/JuMPeR.jl

"""
    setup_set(UncSet, RobustModel, idxs, scens_requested, other_prefs)

Called when a RobustModel is being solved, but before any reformulations or
cuts have been requested. Notifies the `uncset` that it is responsible for
uncertain constraints belonging to the RobustModel with indices `idxs`.
Examples of work that could be done here include transforming the uncertainty
set somehow, or building a cutting plane generating model. Will be called once.
If the user has requested that `Scenario`s be generated at optimality, then
`scens_requested` will be `true` - the uncertainty set may want to generate
the cutting plane model in anticipation of this, even if cutting planes are
not going to be used for solving the problem. `other_prefs` is a dictionary
of keyword arguments passed via the `solve(RobustModel)` function.
"""
setup_set(us::AbstractUncertaintySet, rm::Model, idxs::Vector{Int},
            scens_requested::Bool, other_prefs) =
    error("$(typeof(us)) has not implemented setup_set")


"""
    generate_reform(UncSet, RobustModel, idxs)

Called immediately before the main solve loop (where cutting planes are
generated, if required). Can add anything it wants to the deterministic part
of the RobustModel, but is generally intended for replacing the uncertain
constraints in `idxs` with deterministic equivalents (and possibly adding new
auxiliary variables).
"""
generate_reform(us::AbstractUncertaintySet, rm::Model, idxs::Vector{Int}) =
    error("$(typeof(us)) hasn't implemented generate_reform")


"""
    generate_cut(UncSet, RobustModel, idxs)

Called in the main loop every iteration (continuous variables) or every time
an integer solution is found (discrete variables). Returns a vector of
deterministic constraints which are added to the problem by main solve loop.
Generally intended for generating deterministic versions of the uncertain
constraints in `idxs` as needed.
"""
generate_cut(us::AbstractUncertaintySet, rm::Model, idxs::Vector{Int}) =
    error("$(typeof(us)) hasn't implemented generate_cut")


"""
    generate_scenario(UncSet, RobustModel, idxs)

If requested by the user, this method will be called at optimality. Returns
a `Union{Scenario, Missing}` for each constraint, where that `Scenario` corresponds
to the values of the uncertain parameters that reduce slack in the constraint
the most. If there are multiple such sets of values, the uncertainty set can
select arbitrarily, and if the set cannot provide a scenario it should return
an empty `Union{Scenario, Missing}`.
"""
generate_scenario(us::AbstractUncertaintySet, rm::Model, idxs::Vector{Int}) =
    error("$(typeof(us)) hasn't implemented generate_scenario")


# @constraint support methods for AbstractUncertaintySets
# These do not need to be implemented for all uncertainty sets - only if
# they make sense for the set.
function JuMP.addconstraint(us::AbstractUncertaintySet, c::UncSetConstraint)
    error("$(typeof(us)) hasn't implemented adding constraints on uncertain parameters.")
end
function JuMP.addconstraint(m::JuMP.AbstractModel, c::Array{UncSetConstraint})
    error("The operators <=, >=, and == can only be used to specify scalar constraints. If you are trying to add a vectorized constraint, use the element-wise dot comparison operators (.<=, .>=, or .==) instead")
end
function JuMP.addVectorizedConstraint(m::JuMP.AbstractModel, v::Array{UncSetConstraint})
    map(c->JuMP.addconstraint(m,c), v)
end
function JuMP.addconstraint(us::AbstractUncertaintySet, c::UncSetNormConstraint)
    error("$(typeof(us)) hasn't implemented adding constraints on uncertain parameters.")
end
# Sometimes JuMP[eR] can produce UncConstraints that have no variables - these
# are actually UncSetConstraints, but just haven't been recognized as such.
# To avoid the need for all uncertainty sets to be aware of this, we also
# provide fall backs to detect these.
function JuMP.addconstraint(us::AbstractUncertaintySet, c::UncConstraint)
    if length(c.terms.vars) == 0
        # Pure uncertain constraint
        return JuMP.addconstraint(us, UncSetConstraint(c.terms.constant, c.lb, c.ub))
    end
    # Error, has variables!
    error("Can't add a constraint with decision variables to an uncertainty set!")
end
function JuMP.addVectorizedConstraint(us::AbstractUncertaintySet, v::Array{UncConstraint})
    map(c->JuMP.addconstraint(us,c), v)
end


# Utility functions for common UncertaintySet operations
include("uncsets_util.jl")

# The default BasicUncertaintySet, which handles explicitly provided sets.
include("uncsets_basic.jl")

# BudgetUncertaintySet, based on the set from the "Price of Robustness"
# paper by Bertsimas and Sim.
include("uncsets_budget.jl")
