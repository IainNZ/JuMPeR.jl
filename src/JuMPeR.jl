#-----------------------------------------------------------------------
# JuMPeR  --  JuMP Extension for Robust Optimization
# http://github.com/IainNZ/JuMPeR.jl
#-----------------------------------------------------------------------
# Copyright (c) 2016: Iain Dunning
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
#-----------------------------------------------------------------------
# src/JuMPeR.jl
# Defines the module. Most functionality is `include`d, but this file
# includes RobustModel, UncVarExpr, and UncConstraint.
#-----------------------------------------------------------------------

# __precompile__()

module JuMPeR

importall Base.Operators
import MathProgBase

# Import everything we need from JuMP, so we can build on it
importall JuMP
import JuMP: addconstraint
import JuMP: JuMPConstraint, sense, rhs
import JuMP: GenericAffExpr, GenericRangeConstraint
import JuMP: GenericNorm, GenericNormExpr
import JuMP: IndexedVector, addelt!
import JuMP: JuMPContainer, JuMPDict, JuMPArray

# JuMPeRs exported interface
export RobustModel, @uncertain, uncertainty_set!, getScenario, unc_value


# Forward definition of the abstract parent type for uncertainty sets
# Interface is defined in uncsets.jl
"""
    AbstractUncertaintySet

All uncertainty sets implement the interface defined by AbstractUncertaintySet.
Parent type is JuMP.AbstractModel, to enable JuMP's `@constraint`, etc.
"""
abstract AbstractUncertaintySet <: JuMP.AbstractModel


"""
    RobustModelExt

JuMP extension structure for RobustModel. A RobustModel is a collection of
uncertain parameters (which have similar properties to decision variables),
adaptive variables, and constraints that mix adaptive variables, uncertain
parameters, and normal decision variables. It has a default uncertainty set
that is used for all constraints, unless overriden with a per-constraint set.

Fields:
  Uncertain parameters:
    num_uncs                Number of uncertain parameters in the model
    unc_names               Names
    unc_lower, unc_upper    Lower and upper bounds (across uncertainty sets)
    unc_cat                 Category (normally :Cont, but can be :Int or :Bin)
  Adaptive variables:
    num_adps                Number of adaptive variables in the model
    adp_names               Names
    adp_lower, adp_upper    Lower and upper bounds
    adp_cat                 Category (e.g., :Cont, :Bin)
    adp_policy              Adaptive policy, currently :Static or :Affine
    adp_arguments           Uncertain parameters this variable is a function of
  Constraints:
    unc_constraints         Normal/adapt. variables w/ uncertain parameters
    adapt_constraints       Adaptive variables w/ certain parameters
  Solvers:
    cutsolver               Solver suggested to sets for use when solving
                            cutting plane problems
  Uncertainty sets:
    default_uncset          Default AbstractUncertaintySet used for constraints
    constraint_uncsets      Per-constraint uncertainty sets (or nothing)
  Scenarios:
    scenarios               Stores a Scenario per constraint, if requested.
  Misc:
    solved                  Flags if solved already (to prevent resolves)
"""
type RobustModelExt{S,T,U}
    # Uncertain parameters
    num_uncs::Int
    unc_names::Vector{UTF8String}
    unc_lower::Vector{Float64}
    unc_upper::Vector{Float64}
    unc_cat::Vector{Symbol}
    # Adaptive variables
    num_adps::Int
    adp_names::Vector{UTF8String}
    adp_lower::Vector{Float64}
    adp_upper::Vector{Float64}
    adp_cat::Vector{Symbol}
    adp_policy::Vector{Symbol}
    adp_arguments::Vector{Any}
    # Constraints
    unc_constraints::Vector{S}      # UncConstraint
    adapt_constraints::Vector{T}    # AdaptConstraint
    # Can have different solver for cutting planes
    cutsolver::MathProgBase.SolverInterface.AbstractMathProgSolver
    # Uncertainty sets
    default_uncset::AbstractUncertaintySet
    constraint_uncsets::Vector{Any}
    # Pretty printing magic
    dictList::Vector
    uncDict::Dict{Symbol,Any}
    uncData::ObjectIdDict
    # Scenarios
    scenarios::Vector{Nullable{U}}
    # Misc
    solved::Bool
end
RobustModelExt(cutsolver) =
    RobustModelExt{UncConstraint, AdaptConstraint, Scenario}(
    # Uncertain parameters
    0, UTF8String[],            # num_uncs, unc_names
    Float64[], Float64[],       # unc_lower, unc_upper
    Symbol[],                   # unc_cat
    # Adaptive variables
    0, UTF8String[],            # num_adps, adp_names
    Float64[], Float64[],       # adp_lower, adp_upper
    Symbol[], Symbol[], Any[],  # adp_cat, adp_policy,adp_arguments
    # Constraints
    UncConstraint[],            # unc_constraints
    AdaptConstraint[],          # adapt_constraints
    # cutsolver
    cutsolver,
    # Uncertainty sets
    BasicUncertaintySet(),      # default_uncset
    Any[],                      # constraint_uncsets
    # Pretty printing magic
    Any[],                      # dictList
    Dict{Symbol,Any}(),         # uncDict
    ObjectIdDict(),             # uncData
    # Scenarios
    Nullable{Scenario}[],       # scenarios
    # Misc
    false)                      # solved


"""
    RobustModel(;solver=..., cutsolver=...)

Constructs a extended JuMP model that is the description of an uncertainty set.
The only operations it supports are adding constraints on uncertain parameters
that belong to the parent model (which must be a RobustModel).
"""
function RobustModel(; solver=JuMP.UnsetSolver(),
                    cutsolver=JuMP.UnsetSolver())
    # Create the underlying JuMP model
    m = Model(solver=solver)
    # Add the robust extensions
    m.ext[:JuMPeR] = RobustModelExt(cutsolver)
    # Override the default printing and solving calls
    JuMP.setprinthook(m, print_robust)
    JuMP.setsolvehook(m, solve_robust)
    return m
end


"""
    get_robust(RobustModel)

Internal helper function that simultaneously validates that a JuMP Model is
a JuMPeR RobustModel, and returns the extension object.
"""
function get_robust(m::Model)
    if haskey(m.ext, :JuMPeR)
        return m.ext[:JuMPeR]
    end
    error("This functionality is only available for a JuMPeR RobustModel.")
end


"""
    uncertainty_set!(RobustModel, UncertaintySet)

Sets the default uncertainty set for the model. This set will be used by all
constraints that don't have have a seperately defined uncertainty set.
"""
function uncertainty_set!(m::Model, us::AbstractUncertaintySet)
    get_robust(m).default_uncset = us
    return us
end


# Uncertain, UncExpr, UncSetConstraint, UncSetNorm, UncSetNormConstraint
include("uncertain.jl")


# Adaptive, AdaptExpr, AdaptConstraint, and the type alias JuMPeRVar
include("adaptive.jl")


"""
    UncVarExpr

`∑ⱼ (∑ᵢ aᵢⱼ uᵢ) xⱼ`  --  affine expression of unc. parameters and variables.
"""
typealias UncVarExpr GenericAffExpr{UncExpr,JuMPeRVar}
UncVarExpr() = zero(UncVarExpr)
Base.convert(::Type{UncVarExpr}, c::Number) =
    UncVarExpr(JuMPeRVar[], UncExpr[], UncExpr(c))
Base.convert(::Type{UncVarExpr}, x::JuMPeRVar) =
    UncVarExpr(JuMPeRVar[x],UncExpr[UncExpr(1)], UncExpr())
Base.convert(::Type{UncVarExpr}, aff::AffExpr) =
    UncVarExpr(copy(aff.vars), map(UncExpr,aff.coeffs), UncExpr(aff.constant))
Base.convert(::Type{UncVarExpr}, uaff::UncExpr) =
    UncVarExpr(JuMPeRVar[], UncExpr[], uaff)
function Base.push!(faff::UncVarExpr, new_coeff::Union{Real,Uncertain}, new_var::JuMPeRVar)
    push!(faff.vars, new_var)
    push!(faff.coeffs, UncExpr(new_coeff))
end


"""
    UncConstraint

A constraint with uncertain parameters and variables (i.e., `UncVarExpr`).
"""
typealias UncConstraint GenericRangeConstraint{UncVarExpr}
function JuMP.addconstraint(m::Model, c::UncConstraint; uncset=nothing)
    # Handle the odd special case where there are actually no variables in
    # the constraint - arises from use of macros
    if length(c.terms.vars) == 0
        # Pure uncertain constraint
        return addconstraint(m, UncSetConstraint(c.terms.constant, c.lb, c.ub))
    end
    # Just a regular old constraint
    rmext = get_robust(m)::RobustModelExt
    push!(rmext.unc_constraints, c)
    push!(rmext.constraint_uncsets, uncset)
    return ConstraintRef{Model,UncConstraint}(m, length(rmext.unc_constraints))
end
JuMP.addconstraint(m::Model, c::Array{UncConstraint}) =
    error("The operators <=, >=, and == can only be used to specify scalar constraints. If you are trying to add a vectorized constraint, use the element-wise dot comparison operators (.<=, .>=, or .==) instead")
function JuMP.addVectorizedConstraint(m::Model, v::Array{UncConstraint})
    map(c->addconstraint(m,c), v)
end


"""
    Scenario

A realization of some or all of the uncertain parameters in a model.
"""
type Scenario
    values::Vector{Float64}  # Using NaN as undefined
end


"""
    unc_value(Scenario, Uncertain)

Returns the value of a particular uncertain parameter in the given Scenario.
"""
unc_value(scen::Scenario, u::Uncertain) = scen.values[u.id]


"""
    getScenario(ConstraintRef{RobustModel,UncConstraint})

Get the Scenario for a constraint (as a `Nullable{Scenario}`)
"""
getScenario(uc::ConstraintRef{Model,UncConstraint}) = get_robust(uc.m).scenarios[uc.idx]


# Operator overloads for JuMPeR types
include("operators.jl")

# Adaptive optimization - @adaptive, etc.
include("adaptive/macro.jl")
include("adaptive/expand.jl")

# Define the uncertainty set interface, as well as the BasicUncertaintySet
# and the BudgetUncertaintySet
include("uncsets.jl")

# All functions related to actually solving the model
include("solve.jl")

# Macros for more efficient generation
include("robustmacro.jl")

# Pretty printing
include("print.jl")
include("adaptive/print.jl")

end  # module
