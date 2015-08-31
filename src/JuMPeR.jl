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

isdefined(Base, :__precompile__) && __precompile__()

module JuMPeR

using Compat

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
export UAffExpr, FullAffExpr
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
    0.0)            # solve_time
function RobustModel(; solver=JuMP.UnsetSolver(),
                    cutsolver=JuMP.UnsetSolver())
    # Create the underlying JuMP model
    m = Model(solver=solver)
    # Add the robust extensions
    m.ext[:JuMPeR] = RobustData(cutsolver)
    # Override the default printing and solving calls
    JuMP.setPrintHook(m, _print_robust)
    JuMP.setSolveHook(m, _solve_robust)
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
function Uncertain(m::Model, lower::Number, upper::Number, cat::Symbol, name::String)
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
Base.zero(::Type{Uncertain}) = UAffExpr()
Base.zero(::Uncertain) = zero(Uncertain)
Base.one(::Type{Uncertain}) = UAffExpr(1)
Base.one(::Uncertain) = one(Uncertain)
Base.isequal(u1::Uncertain, u2::Uncertain) = (u1.m === u2.m) && isequal(u1.id, u2.id)
getNumUncs(m::Model) = getRobust(m).numUncs


#-----------------------------------------------------------------------
# UAffExpr   ∑ᵢ aᵢ uᵢ
typealias UAffExpr GenericAffExpr{Float64,Uncertain}

UAffExpr() = zero(UAffExpr)
UAffExpr(x::Union(Number,Uncertain)) = convert(UAffExpr, x)
UAffExpr(c::Number,u::Uncertain) = UAffExpr(Uncertain[u],Float64[c],0.0)
Base.convert(::Type{UAffExpr}, u::Uncertain) = UAffExpr(Uncertain[u],Float64[1],0.0)
Base.convert(::Type{UAffExpr}, c::Number)    = UAffExpr(Uncertain[ ],Float64[ ],  c)
# aff_to_uaff
# Useful for oracles. Given a UAffExpr and a list of variables, create an
# AffExpr such that Uncertain(i) maps to Variable(i), where i is the index,
# in the new expression.
uaff_to_aff(uaff::UAffExpr, x::Vector{Variable}) =
    AffExpr(Variable[x[up.id] for up in uaff.vars],
            copy(uaff.coeffs), uaff.constant)


#-----------------------------------------------------------------------
# FullAffExpr   ∑ⱼ (∑ᵢ aᵢⱼ uᵢ) xⱼ
typealias FullAffExpr GenericAffExpr{UAffExpr,Variable}

FullAffExpr() = zero(FullAffExpr)
Base.convert(::Type{FullAffExpr}, c::Number) =
    FullAffExpr(Variable[], UAffExpr[], UAffExpr(c))
Base.convert(::Type{FullAffExpr}, x::Variable) =
    FullAffExpr(Variable[x],UAffExpr[UAffExpr(1)], UAffExpr())
Base.convert(::Type{FullAffExpr}, aff::AffExpr) =
    FullAffExpr(copy(aff.vars), map(UAffExpr,aff.coeffs), UAffExpr(aff.constant))
Base.convert(::Type{FullAffExpr}, uaff::UAffExpr) =
    FullAffExpr(Variable[], UAffExpr[], uaff)

function Base.push!(faff::FullAffExpr, new_coeff::Union(Real,Uncertain), new_var::Variable)
    push!(faff.vars, new_var)
    push!(faff.coeffs, UAffExpr(new_coeff))
end


#-----------------------------------------------------------------------
# UncSetConstraint      A constraint with just uncertain parameters
typealias UncSetConstraint GenericRangeConstraint{UAffExpr}
addConstraint(m::Model, c::UncSetConstraint) = push!(getRobust(m).uncertaintyset, c)
addConstraint(m::Model, c::Array{UncSetConstraint}) =
    error("The operators <=, >=, and == can only be used to specify scalar constraints. If you are trying to add a vectorized constraint, use the element-wise dot comparison operators (.<=, .>=, or .==) instead")
function addVectorizedConstraint(m::Model, v::Array{UncSetConstraint})
    map(c->addConstraint(m,c), v)
    v
end


#-----------------------------------------------------------------------
# UncConstraint         A constraint with variables and uncertains
typealias UncConstraint GenericRangeConstraint{FullAffExpr}
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
JuMP._build_norm(Lp, terms::Vector{UAffExpr}) = UncSetNorm{Lp}(terms)

type UncNormConstraint{P} <: JuMPConstraint
    normexpr::GenericNormExpr{P,Float64,Uncertain}
end

function addConstraint(m::Model, c::UncNormConstraint)
    push!(getRobust(m).normconstraints,c)
    getRobust(m).normconstraints[end]
end

addEllipseConstraint(m::Model, vec::Vector, Gamma::Real) =
    error("""addEllipseConstraint not supported as of JuMPeR v0.2.
             Please use, e.g., @addConstraint(m, norm(x) <= Γ)
                               @addConstraint(m, norm2{x[i],i=1:n} <= Γ""")


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