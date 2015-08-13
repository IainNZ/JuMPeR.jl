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

module JuMPeR

using Compat

importall Base.Operators

# Import everything we need from JuMP, so we can build on it
importall JuMP
import JuMP: GenericAffExpr, GenericRangeConstraint
import JuMP: GenericNorm, GenericNormExpr
import JuMP: JuMPConstraint
import JuMP: sense, rhs
import JuMP: addVectorizedConstraint
import JuMP: IndexedVector, addelt! #, isexpr
import JuMP: JuMPContainer, JuMPDict, JuMPArray

# JuMPeRs exported interface
export RobustModel, getNumUncs
export setDefaultOracle!
export Uncertain, @defUnc, addEllipseConstraint
export UAffExpr, FullAffExpr
export UncConstraint, UncSetConstraint, EllipseConstraint
# Deprecated
export printRobust, solveRobust



#############################################################################
# RobustData contains all extensions to the base JuMP model type
type RobustData
    # Variable-Uncertain mixed constraints
    uncertainconstr
    # Oracles associated with each uncertainconstr
    oracles
    # Uncertain-only constraints
    uncertaintyset
    normconstraints
    
    # Uncertainty data
    numUncs::Int
    uncNames::Vector{String}
    uncLower::Vector{Float64}
    uncUpper::Vector{Float64}
    uncCat::Vector{Symbol}

    # Adaptability
    adapt_type::Dict{Int,Symbol}
    adapt_on::Dict{Int,Vector}

    defaultOracle

    # Active cuts
    activecuts

    # Can have different solver for cutting planes
    cutsolver

    # For pretty printing
    dictList::Vector

    # Provided scenarios
    scenarios::Vector

    solve_time::Float64
end
RobustData(cutsolver) = RobustData(Any[],Any[],Any[],Any[],
                            0,String[],Float64[],Float64[],Symbol[],
                            Dict{Int,Symbol}(), Dict{Int,Vector}(),
                            GeneralOracle(), Any[],
                            cutsolver,JuMP.JuMPContainer[],Any[], 0.0)

function RobustModel(;solver=JuMP.UnsetSolver(),cutsolver=JuMP.UnsetSolver())
    m = Model(solver=solver)
    m.ext[:Robust] = RobustData(cutsolver)
    JuMP.setPrintHook(m, _print_robust)
    JuMP.setSolveHook(m, _solve_robust)
    return m
end

function getRobust(m::Model)
    if haskey(m.ext, :Robust)
        return m.ext[:Robust]
    else
        error("This functionality is only available for RobustModels")
    end
end

getNumUncs(m::Model) = getRobust(m).numUncs

#############################################################################
# Uncertain
# Similar to JuMP.Variable, has an reference back to the model and an id num
type Uncertain <: JuMP.AbstractJuMPScalar
    m::Model
    unc::Int
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
getLower(u::Uncertain) = getRobust(u.m).uncLower[u.unc]
getUpper(u::Uncertain) = getRobust(u.m).uncUpper[u.unc]
getName(u::Uncertain) = unc_str(REPLMode, u.m, u.unc)
getCategory(u::Uncertain) = getRobust(u.m).uncCat[u.unc]
Base.zero(::Type{Uncertain}) = UAffExpr()
Base.zero(::Uncertain) = zero(Uncertain)
Base.one(::Type{Uncertain}) = UAffExpr(1)
Base.one(::Uncertain) = one(Uncertain)
Base.isequal(u1::Uncertain, u2::Uncertain) = isequal(u1.unc, u2.unc)
# Generic matrix operators in Base expect to be able to find
# a type to use for the result, but that type needs to have a zero
# so we can't leave it as Any
#Base.promote_rule{T<:Real}(::Type{Uncertain},::Type{T}) = UAffExpr

#############################################################################
# Uncertain Affine Expression class
typealias UAffExpr GenericAffExpr{Float64,Uncertain}

UAffExpr() = UAffExpr(Uncertain[],Float64[],0.)
UAffExpr(c::Real) = UAffExpr(Uncertain[],Float64[],float(c))
UAffExpr(u::Uncertain) = UAffExpr([u],[1.],0.)
UAffExpr(c::Real,u::Uncertain,constant=0) = UAffExpr([u],[float(c)],constant)
Base.convert(::Type{UAffExpr}, u::Uncertain) = UAffExpr(u)
Base.convert(::Type{UAffExpr}, c::Number) = UAffExpr(c)
# aff_to_uaff
# Useful for oracles. Given a UAffExpr and a list of variables, create an
# AffExpr such that Uncertain(i) maps to Variable(i), where i is the index,
# in the new expression.
uaff_to_aff(uaff::UAffExpr, x::Vector{Variable}) =
    AffExpr(Variable[x[up.unc] for up in uaff.vars],
            copy(uaff.coeffs), uaff.constant)

#############################################################################
# Full Affine Expression class
# Todo: better name. In my other robust modelling tools I called it
# something like this, but the catch then was that there we only two types of
# affexpr - the one with UAffExpr coefficients = Full, and the UAffExpr itself
typealias FullAffExpr GenericAffExpr{UAffExpr,Variable}

FullAffExpr() = FullAffExpr(Variable[], UAffExpr[], UAffExpr())
Base.zero(a::Type{FullAffExpr}) = FullAffExpr()
Base.zero(a::FullAffExpr) = zero(typeof(a))
Base.convert(::Type{FullAffExpr}, x::Variable) = FullAffExpr([x],[UAffExpr(1)], UAffExpr())
Base.convert(::Type{FullAffExpr}, aff::AffExpr) =
    FullAffExpr(aff.vars,map(UAffExpr,aff.coeffs), UAffExpr(aff.constant))
function Base.push!(faff::FullAffExpr, new_coeff::UAffExpr, new_var::Variable)
    push!(faff.vars, new_var)
    push!(faff.coeffs, new_coeff)
end
function Base.push!(faff::FullAffExpr, new_coeff::Union(Real,Uncertain), new_var::Variable)
    push!(faff.vars, new_var)
    push!(faff.coeffs, UAffExpr(new_coeff))
end

#############################################################################
# UncSetConstraint      Just uncertainties
typealias UncSetConstraint GenericRangeConstraint{UAffExpr}
addConstraint(m::Model, c::UncSetConstraint) = push!(getRobust(m).uncertaintyset, c)
addConstraint(m::Model, c::Array{UncSetConstraint}) =
    error("The operators <=, >=, and == can only be used to specify scalar constraints. If you are trying to add a vectorized constraint, use the element-wise dot comparison operators (.<=, .>=, or .==) instead")
addVectorizedConstraint(m::Model, v::Array{UncSetConstraint}) = map(c->addConstraint(m,c), v)

# UncConstraint         Mix of variables and uncertains
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
addVectorizedConstraint(m::Model, v::Array{UncConstraint}) = map(c->addConstraint(m,c), v)

#############################################################################
# Norms
import JuMP: _build_norm
typealias UncSetNorm{Typ} GenericNorm{Typ,Float64,Uncertain}
_build_norm(Lp, terms::Vector{UAffExpr}) = UncSetNorm{Lp}(terms)

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

#############################################################################
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

#############################################################################
end  # module
#############################################################################