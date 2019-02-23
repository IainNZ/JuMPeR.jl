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

import MathProgBase
using JuMP  # So we can build on it, but prefer explicit qualification
import JuMP: JuMPContainer, GenericAffExpr, GenericNorm, GenericNormExpr, getname
import LinearAlgebra: dot, norm

# JuMPeRs exported interface
export RobustModel, @uncertain, @adaptive,
        setuncertaintyset, getscenario, uncvalue


# Forward definition of the abstract parent type for uncertainty sets
# Interface is defined in uncsets.jl
"""
    AbstractUncertaintySet

All uncertainty sets implement the interface defined by AbstractUncertaintySet.
Parent type is JuMP.AbstractModel, to enable JuMP's `@constraint`, etc.
"""
abstract type AbstractUncertaintySet <: JuMP.AbstractModel end


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
mutable struct RobustModelExt{S,T,U}
    # Uncertain parameters
    num_uncs::Int
    unc_names::Vector{String}
    unc_lower::Vector{Float64}
    unc_upper::Vector{Float64}
    unc_cat::Vector{Symbol}
    # Adaptive variables
    num_adps::Int
    adp_names::Vector{String}
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
    # Scenarios
    scenarios::Vector{Union{U, Missing}}
    # Misc
    solved::Bool
    # Pretty printing magic
    dictList::Vector
    uncDict::Dict{Symbol,Any}
    uncData::IdDict{Any, Any}
end
RobustModelExt(cutsolver) =
    RobustModelExt{UncConstraint, AdaptConstraint, Scenario}(
    # Uncertain parameters
    0, String[],                # num_uncs, unc_names
    Float64[], Float64[],       # unc_lower, unc_upper
    Symbol[],                   # unc_cat
    # Adaptive variables
    0, String[],                # num_adps, adp_names
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
    # Scenarios
    Union{Scenario, Missing}[], # scenarios
    # Misc
    false,                      # solved
    # Pretty printing magic
    Any[],                      # dictList
    Dict{Symbol,Any}(),         # uncDict
    IdDict{Any, Any}())         # uncData


"""
    RobustModel(;solver=..., cutsolver=...)

Constructs a JuMP model with a RobustModelExt extension, and new solve and
print functions via the respective hooks.
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
    setuncertaintyset(RobustModel, UncertaintySet)

Sets the default uncertainty set for the model. This set will be used by all
constraints that don't have have a seperately defined uncertainty set.
"""
function setuncertaintyset(m::Model, us::AbstractUncertaintySet)
    get_robust(m).default_uncset = us
    return us
end

# TODO: properly register uncset constraints
JuMP.registercon(m::AbstractUncertaintySet, conname, value) = value


# Uncertain, UncExpr, UncSetConstraint, UncSetNormConstraint
include("uncertain.jl")


# Adaptive, AdaptExpr, AdaptConstraint, and the type alias JuMPeRVar
include("adaptive.jl")


"""
    UncVarExpr

`∑ⱼ (∑ᵢ aᵢⱼ uᵢ) xⱼ`  --  affine expression of unc. parameters and variables.
"""
UncVarExpr = JuMP.GenericAffExpr{UncExpr,JuMPeRVar}
Base.convert(::Type{UncVarExpr}, c::Number) =
    UncVarExpr(JuMPeRVar[], UncExpr[], UncExpr(c))
Base.convert(::Type{UncVarExpr}, x::JuMPeRVar) =
    UncVarExpr(JuMPeRVar[x],UncExpr[UncExpr(1)], UncExpr())
Base.convert(::Type{UncVarExpr}, aff::AffExpr) =
    UncVarExpr(copy(aff.vars), convert.(UncExpr,aff.coeffs), UncExpr(aff.constant))
Base.convert(::Type{UncVarExpr}, uaff::UncExpr) =
    UncVarExpr(JuMPeRVar[], UncExpr[], uaff)
JuMP.GenericAffExpr{U,V}() where {U<:UncExpr,V<:JuMPeRVar} = zero(UncVarExpr)
JuMP.GenericAffExpr{U,V}(x::Union{JuMPeRVar,AffExpr,UncExpr}) where {U<:UncExpr,V<:JuMPeRVar} = convert(UncExpr, x)
function Base.push!(faff::UncVarExpr, new_coeff::Union{Real,Uncertain}, new_var::JuMPeRVar)
    push!(faff.vars, new_var)
    push!(faff.coeffs, UncExpr(new_coeff))
end


"""
    UncConstraint

A constraint with uncertain parameters and variables (i.e., `UncVarExpr`).
"""
UncConstraint = JuMP.GenericRangeConstraint{UncVarExpr}
function JuMP.addconstraint(m::Model, c::UncConstraint; uncset=nothing)
    # Handle the odd special case where there are actually no variables in
    # the constraint - arises from use of macros
    if length(c.terms.vars) == 0
        # Pure uncertain constraint
        return JuMP.addconstraint(m, UncSetConstraint(c.terms.constant, c.lb, c.ub))
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
    map(c->JuMP.addconstraint(m,c), v)
end


"""
    Scenario

A realization of some or all of the uncertain parameters in a model.
"""
mutable struct Scenario
    values::Vector{Float64}  # Using NaN as undefined
end


"""
    uncvalue(Scenario, Uncertain)

Returns the value of a particular uncertain parameter in the given Scenario.
"""
uncvalue(scen::Scenario, u::Uncertain) = scen.values[u.id]


"""
    getscenario(ConstraintRef{RobustModel,UncConstraint})

Get the Scenario for a constraint (as a `Union{Scenario, Missing}`)
"""
getscenario(uc::ConstraintRef{Model,UncConstraint}) = get_robust(uc.m).scenarios[uc.idx]


# Operator overloads for JuMPeR types
include("operators.jl")

# Adaptive optimization support
include("expand.jl")

# Define the uncertainty set interface, as well as the BasicUncertaintySet
# and the BudgetUncertaintySet
include("uncsets.jl")

# All functions related to actually solving the model
include("solve.jl")

# Macros for more efficient generation
include("robustmacro.jl")

# Pretty printing
include("print.jl")

end  # module
