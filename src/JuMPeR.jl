#-----------------------------------------------------------------------
# JuMPeR  --  JuMP Extension for Robust Optimization
# http://github.com/IainNZ/JuMPeR.jl
#-----------------------------------------------------------------------
# Copyright (c) 2015: Iain Dunning
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
#-----------------------------------------------------------------------

module JuMPeR

using Compat

importall Base.Operators

# Import everything we need from JuMP, so we can build on it
importall JuMP
import JuMP: GenericAffExpr, GenericRangeConstraint
import JuMP: sense, rhs
import JuMP.IndexedVector, JuMP.addelt!, JuMP.isexpr
import JuMP: JuMPContainer, JuMPDict, JuMPArray
import JuMP.@gendict
import JuMP: assert_validmodel, validmodel, esc_nonconstant
import JuMP: getloopedcode, buildrefsets, getname
import JuMP: addVectorizedConstraint

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
getName(u::Uncertain) = unc_str(REPLMode, u.m, u.unc)
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
UAffExpr(u::Uncertain, c::Real) = UAffExpr([u],[float(c)],0.)
Base.convert(::Type{UAffExpr}, u::Uncertain) = UAffExpr(u)
Base.convert(::Type{UAffExpr}, c::Number) = UAffExpr(c)

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
# EllipseConstraint
# Capture uncertainty set constraints of the form  || F u + g ||_2 <= Gamma
type EllipseConstraint <: JuMP.JuMPConstraint
    m::Model
    F::Array{Float64, 2}
    u::Array{Int, 1}
    g::Array{Float64, 1}
    Gamma::Float64
end

# build_ellipse_constraint
# Given || vec ||_2 <= Gamma, return an EllipseConstraint by expanding 
# `vec` out to its full `Fu+g` form.
function build_ellipse_constraint(m::Model, vec::Vector, Gamma::Float64)
    
    # In the first pass we determine a unique set of uncertainties
    # present so we can allocate the correct size F and u
    unc_map  = Dict{Int,Int}()
    rev_map  = Dict{Int,Int}()
    function record(v::Uncertain)
        if !(v.unc in keys(unc_map))
            unc_map[v.unc] = length(unc_map) + 1
            rev_map[length(unc_map)] = v.unc
        end
    end
    record(v::UAffExpr) = map(record, v.vars)
    record(v) = ArgumentError()
    map(record, vec)

    # Create F and g
    num_uncs  = length(unc_map)
    num_terms = length(vec)
    F = zeros(num_terms, num_uncs)
    g = zeros(num_terms)
    for (i,v) in enumerate(vec)
        if typeof(v) <: UAffExpr
            g[i] = v.constant
            for (j,u) in enumerate(v.vars)
                F[i,unc_map[u.unc]] += v.coeffs[j]
            end
        elseif typeof(v) <: Uncertain
            F[i,unc_map[v.unc]] += 1.0
        end
    end

    return EllipseConstraint(m,F,[rev_map[i] for i=1:num_uncs],g,Gamma)
end

function addEllipseConstraint(m::Model, vec::Vector, Gamma::Real)
    push!(getRobust(m).normconstraints, build_ellipse_constraint(m,vec,float(Gamma)))
    getRobust(m).normconstraints[end]
end

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