#-----------------------------------------------------------------------
# JuMPeR  --  JuMP Extension for Robust Optimization
# http://github.com/IainNZ/JuMPeR.jl
#-----------------------------------------------------------------------
# Copyright (c) 2015: Iain Dunning
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
#-----------------------------------------------------------------------
# src/JuMPeR.jl
# Defines the module and main types: RobustData, Uncertain, and the
# various constraints and norms.
#-----------------------------------------------------------------------

__precompile__()

module JuMPeR

importall Base.Operators

# Import everything we need from JuMP, so we can build on it
importall JuMP
import JuMP: JuMPConstraint, sense, rhs
import JuMP: GenericAffExpr, GenericRangeConstraint
import JuMP: GenericNorm, GenericNormExpr
import JuMP: addVectorizedConstraint
import JuMP: IndexedVector, addelt!
import JuMP: JuMPContainer, JuMPDict, JuMPArray

# JuMPeRs exported interface
export RobustModel, getNumUncs
export setDefaultOracle!
export Uncertain, @defUnc, addEllipseConstraint
export UncExpr, UncAffExpr
export UncConstraint, UncSetConstraint, EllipseConstraint


#-----------------------------------------------------------------------
# RobustData contains all extensions to the base JuMP model type
type RobustData
    # Variable-Uncertain mixed constraints
    uncertainconstr::Vector
    # Oracles associated with each uncertainconstr
    oracles::Vector
    # Uncertain-only constraints
    uncertaintyset::Vector
    normconstraints::Vector

    # Uncertainty data
    numUncs::Int
    uncNames::Vector{UTF8String}
    uncLower::Vector{Float64}
    uncUpper::Vector{Float64}
    uncCat::Vector{Symbol}
    defaultOracle

    # Active cuts
    activecuts

    # Can have different solver for cutting planes
    cutsolver

    # For pretty printing
    dictList::Vector
    uncDict::Dict{Symbol,Any}
    uncData::ObjectIdDict

    # Provided scenarios
    scenarios::Vector

    solve_time::Float64

    # Custom lazy constraint callbacks, which need special handling due
    # to normal lazy constraint callbacks attaching to the outer RobustModel
    lazy_callbacks::Vector{Any}
end
RobustData(cutsolver) = RobustData(
    Any[],          # uncertainconstr
    Any[],          # oracles
    Any[],          # uncertaintyset
    Any[],          # normconstraints
    0,              # numUncs
    UTF8String[],   # uncNames
    Float64[],      # uncLower
    Float64[],      # uncUpper
    Symbol[],       # uncCat
    GeneralOracle(),# defaultOracle
    Any[],          # activecuts
    cutsolver,      # cutsolver
    Any[],          # dictList
    Dict{Symbol,Any}(), # uncDict
    ObjectIdDict(), # uncData
    Any[],          # scenarios
    0.0,            # solve_time
    Any[])          # lazy_callbacks
function RobustModel(; solver=JuMP.UnsetSolver(),
                    cutsolver=JuMP.UnsetSolver())
    # Create the underlying JuMP model
    m = Model(solver=solver)
    # Add the robust extensions
    m.ext[:JuMPeR] = RobustData(cutsolver)
    # Override the default printing and solving calls
    JuMP.setPrintHook(m, print_robust)
    JuMP.setSolveHook(m, solve_robust)
    return m
end

function getRobust(m::Model)
    haskey(m.ext, :JuMPeR) && return m.ext[:JuMPeR]
    error("This functionality is only available for JuMPeR RobustModels.")
end


#-----------------------------------------------------------------------
# Uncertain
# Similar to JuMP.Variable, has an reference back to the model and an id num
type Uncertain <: JuMP.AbstractJuMPScalar
    m::Model
    id::Int
end
function Uncertain(m::Model, lower::Number, upper::Number, cat::Symbol, name::AbstractString)
    robdata = getRobust(m)
    robdata.numUncs += 1
    push!(robdata.uncNames, name)
    push!(robdata.uncLower, convert(Float64, lower))
    push!(robdata.uncUpper, convert(Float64, upper))
    push!(robdata.uncCat, cat)
    return Uncertain(m, robdata.numUncs)
end
Uncertain(m::Model, lower::Number, upper::Number, cat::Symbol) = Uncertain(m,lower,upper,cat,"")
getLower(u::Uncertain) = getRobust(u.m).uncLower[u.id]
getUpper(u::Uncertain) = getRobust(u.m).uncUpper[u.id]
getName(u::Uncertain) = unc_str(REPLMode, u.m, u.id)
getCategory(u::Uncertain) = getRobust(u.m).uncCat[u.id]
Base.zero(::Type{Uncertain}) = UncExpr()
Base.zero(::Uncertain) = zero(Uncertain)
Base.one(::Type{Uncertain}) = UncExpr(1)
Base.one(::Uncertain) = one(Uncertain)
Base.isequal(u1::Uncertain, u2::Uncertain) = (u1.m === u2.m) && isequal(u1.id, u2.id)
getNumUncs(m::Model) = getRobust(m).numUncs


#-----------------------------------------------------------------------
# UncExpr   ∑ᵢ aᵢ uᵢ
typealias UncExpr GenericAffExpr{Float64,Uncertain}

UncExpr() = zero(UncExpr)
UncExpr(x::Union{Number,Uncertain}) = convert(UncExpr, x)
UncExpr(c::Number,u::Uncertain) = UncExpr(Uncertain[u],Float64[c],0.0)
Base.convert(::Type{UncExpr}, u::Uncertain) = UncExpr(Uncertain[u],Float64[1],0.0)
Base.convert(::Type{UncExpr}, c::Number)    = UncExpr(Uncertain[ ],Float64[ ],  c)
# aff_to_uaff
# Useful for oracles. Given a UncExpr and a list of variables, create an
# AffExpr such that Uncertain(i) maps to Variable(i), where i is the index,
# in the new expression.
uaff_to_aff(uaff::UncExpr, x::Vector{Variable}) =
    AffExpr(Variable[x[up.id] for up in uaff.vars],
            copy(uaff.coeffs), uaff.constant)

# UncSetConstraint      A constraint with just uncertain parameters
typealias UncSetConstraint GenericRangeConstraint{UncExpr}
addConstraint(m::Model, c::UncSetConstraint) = push!(getRobust(m).uncertaintyset, c)
addConstraint(m::Model, c::Array{UncSetConstraint}) =
    error("The operators <=, >=, and == can only be used to specify scalar constraints. If you are trying to add a vectorized constraint, use the element-wise dot comparison operators (.<=, .>=, or .==) instead")
function addVectorizedConstraint(m::Model, v::Array{UncSetConstraint})
    map(c->addConstraint(m,c), v)
    v
end


#-----------------------------------------------------------------------
# UncAffExpr   ∑ⱼ (∑ᵢ aᵢⱼ uᵢ) xⱼ
typealias UncAffExpr GenericAffExpr{UncExpr,Variable}

UncAffExpr() = zero(UncAffExpr)
Base.convert(::Type{UncAffExpr}, c::Number) =
    UncAffExpr(Variable[], UncExpr[], UncExpr(c))
Base.convert(::Type{UncAffExpr}, x::Variable) =
    UncAffExpr(Variable[x],UncExpr[UncExpr(1)], UncExpr())
Base.convert(::Type{UncAffExpr}, aff::AffExpr) =
    UncAffExpr(copy(aff.vars), map(UncExpr,aff.coeffs), UncExpr(aff.constant))
Base.convert(::Type{UncAffExpr}, uaff::UncExpr) =
    UncAffExpr(Variable[], UncExpr[], uaff)

function Base.push!(faff::UncAffExpr, new_coeff::Union{Real,Uncertain}, new_var::Variable)
    push!(faff.vars, new_var)
    push!(faff.coeffs, UncExpr(new_coeff))
end

# UncConstraint         A constraint with variables and uncertains
typealias UncConstraint GenericRangeConstraint{UncAffExpr}
function addConstraint(m::Model, c::UncConstraint, w=nothing)
    rd = getRobust(m)
    # Handle the odd special case where there are actually no variables in
    # the constraint - arises from use of macros
    if length(c.terms.vars) == 0
        # Pure uncertain
        @assert w == nothing
        return addConstraint(m, UncSetConstraint(c.terms.constant, c.lb, c.ub))
    end
    # Just a regular old constraint
    push!(rd.uncertainconstr,c)
    push!(rd.oracles, w)
    push!(rd.activecuts, Any[])
    return ConstraintRef{UncConstraint}(m,length(rd.uncertainconstr))
end
addConstraint(m::Model, c::Array{UncConstraint}) =
    error("The operators <=, >=, and == can only be used to specify scalar constraints. If you are trying to add a vectorized constraint, use the element-wise dot comparison operators (.<=, .>=, or .==) instead")
function addVectorizedConstraint(m::Model, v::Array{UncConstraint})
    map(c->addConstraint(m,c), v)
    v
end


#-----------------------------------------------------------------------
# Norms of uncertain parameters
typealias UncSetNorm{Typ} GenericNorm{Typ,Float64,Uncertain}
JuMP._build_norm(Lp, terms::Vector{UncExpr}) = UncSetNorm{Lp}(terms)

type UncNormConstraint{P} <: JuMPConstraint
    normexpr::GenericNormExpr{P,Float64,Uncertain}
end

function addConstraint(m::Model, c::UncNormConstraint)
    push!(getRobust(m).normconstraints,c)
    getRobust(m).normconstraints[end]
end

#-----------------------------------------------------------------------
# Provide a way to provide a lazy constraint callback for JuMPeR models
# The key problem is that the solution to the internal reformulated
# model is not made available to the JuMPeR model for use in a normal
# JuMP-style callback, and that constraints added to the JuMPeR model
# aren't added to the internal model. Thus, this function provides a
# mechanism to avoid that.
# "Special" lazy callbacks can be written exactly like normal JuMP
# lazy callbacks, except that instead of using @addLazyConstraint, the
# function should return an array of constraints to add. For example,
# function lazycut(cb)
#   x_val = getValue(x)
#   if sum(x_val) > 5
#     return [@LinearConstraint(sum(x) <= 5)]
#   else
#     return nothing
#   end
# end
function addSpecialLazyCallback(m::Model, f::Function)
    push!(getRobust(m).lazy_callbacks, f)
end

#-----------------------------------------------------------------------
# Scenarios
include("scenario.jl")

# Operator overloads
include("robustops.jl")

# All functions related to actual solution
include("solve.jl")

# Oracles... to be name changed
include("oracle.jl")

# Macros for more efficient generation
include("robustmacro.jl")

# Pretty printing
include("print.jl")

# Graph algorithms
include("graph.jl")


#-----------------------------------------------------------------------
end  # module
#-----------------------------------------------------------------------
