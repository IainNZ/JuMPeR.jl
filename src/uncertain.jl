#-----------------------------------------------------------------------
# JuMPeR  --  JuMP Extension for Robust Optimization
# http://github.com/IainNZ/JuMPeR.jl
#-----------------------------------------------------------------------
# Copyright (c) 2016: Iain Dunning
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
#-----------------------------------------------------------------------
# src/uncertain.jl
# Defines the uncertain parameter type Uncertain, and expressions and
# constraints of Uncertains: Uncertain, UncExpr, UncSetConstraint,
# UncSetNorm, UncSetNormConstraint.
# Included by src/JuMPeR.jl
#-----------------------------------------------------------------------


"""
    Uncertain

Uncertain parameter, implemented much the same as a JuMP.Variable. It belongs
to a RobustModel, not to any particular UncertaintySet.
"""
type Uncertain <: JuMP.AbstractJuMPScalar
    m::Model
    id::Int
end
function Uncertain(m::Model, lower::Number, upper::Number, cat::Symbol, name::AbstractString)
    robdata = get_robust(m)::RobustModelExt
    robdata.num_uncs += 1
    push!(robdata.unc_names, convert(UTF8String, name))
    push!(robdata.unc_lower, convert(Float64, lower))
    push!(robdata.unc_upper, convert(Float64, upper))
    push!(robdata.unc_cat, cat)
    return Uncertain(m, robdata.num_uncs)
end
Base.zero(::Type{Uncertain}) = UncExpr()
Base.zero(::Uncertain) = zero(Uncertain)
Base.one(::Type{Uncertain}) = UncExpr(1)
Base.one(::Uncertain) = one(Uncertain)
Base.isequal(u1::Uncertain, u2::Uncertain) = (u1.m === u2.m) && isequal(u1.id, u2.id)


"""
    UncExpr

`∑ᵢ aᵢ uᵢ`  --  affine expression of uncertain parameters and numbers.
"""
typealias UncExpr GenericAffExpr{Float64,Uncertain}
Base.convert(::Type{UncExpr}, u::Uncertain) = UncExpr(Uncertain[u],Float64[1],0.0)
Base.convert(::Type{UncExpr}, c::Number)    = UncExpr(Uncertain[ ],Float64[ ],  c)
UncExpr() = zero(UncExpr)
UncExpr(x::Union{Number,Uncertain}) = convert(UncExpr, x)
UncExpr(c::Number,u::Uncertain) = UncExpr(Uncertain[u],Float64[c],0.0)


"""
    UncSetConstraint

A constraint with just uncertain parameters and numbers (i.e., `UncExpr`).
"""
typealias UncSetConstraint GenericRangeConstraint{UncExpr}
function JuMP.addconstraint(m::Model, c::UncSetConstraint)
    # If m is a RobustModel, we add it to the default uncertainty set
    rmext = get_robust(m)::RobustModelExt
    return addconstraint(rmext.default_uncset, c)
end
function JuMP.addconstraint(m::Model, c::Array{UncSetConstraint})
    error("The operators <=, >=, and == can only be used to specify scalar constraints. If you are trying to add a vectorized constraint, use the element-wise dot comparison operators (.<=, .>=, or .==) instead")
end
function JuMP.addVectorizedConstraint(m::Model, v::Array{UncSetConstraint})
    map(c->addconstraint(m,c), v)
end


"""
    UncSetNorm

A norm on uncertain parameters.
"""
typealias UncSetNorm{Typ} GenericNorm{Typ,Float64,Uncertain}
JuMP._build_norm(Lp, terms::Vector{UncExpr}) = UncSetNorm{Lp}(terms)


"""
    UncSetNormConstraint

A constraint involving a norm of uncertain parameters.
"""
type UncSetNormConstraint{P} <: JuMPConstraint
    normexpr::GenericNormExpr{P,Float64,Uncertain}
end
function JuMP.addconstraint(m::Model, c::UncSetNormConstraint)
    # If m is a RobustModel, we add it to the default uncertainty set
    rmext = get_robust(m)::RobustModelExt
    return addconstraint(rmext.default_uncset, c)
end
