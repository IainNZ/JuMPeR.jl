#-----------------------------------------------------------------------
# JuMPeR  --  JuMP Extension for Robust Optimization
# http://github.com/IainNZ/JuMPeR.jl
#-----------------------------------------------------------------------
# Copyright (c) 2016: Iain Dunning
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
#-----------------------------------------------------------------------
# src/adaptive.jl
# Defines the adaptive variable type Adaptive, and expressions and
# constraints of Adaptive: Adaptive, AdaptExpr, AdaptConstraint, and
# the type alias JuMPeRVar.
# Included by src/JuMPeR.jl
#-----------------------------------------------------------------------

"""
    Adaptive

An adaptive variable, a variable whose value depends on the realized
values of the uncertain parameters.
"""
type Adaptive <: JuMP.AbstractJuMPScalar
    m::Model
    id::Int
end
function Adaptive(m::Model, lower::Real, upper::Real,
                    cat::Symbol, name::AbstractString,
                    policy::Symbol, depends_on::Any)
    rmext = get_robust(m)
    rmext.num_adps += 1
    push!(rmext.adp_lower,  lower)
    push!(rmext.adp_upper,  upper)
    push!(rmext.adp_cat,    cat)
    push!(rmext.adp_names,  name)
    push!(rmext.adp_policy, policy)
    push!(rmext.adp_arguments, depends_on)
    return Adaptive(m, rmext.num_adps)
end
Base.zero(::Type{Adaptive}) = AdaptExpr()
Base.zero(     ::Adaptive)  = zero(Adaptive)
Base.one(::Type{Adaptive})  = AdaptExpr(1)
Base.one(     ::Adaptive)   = one(Adaptive)
Base.isequal(a::Adaptive, b::Adaptive) = (a.m === b.m) && (a.id == b.id)
getname(x::Adaptive) = get_robust(x.m).adp_names[x.id]


"""
    JuMPeRVar

Either a plain JuMP Variable, or a JuMPeR Adaptive variable.
"""
typealias JuMPeRVar Union{Variable,Adaptive}


"""
    AdaptExpr

`∑ᵢ aᵢ vᵢ`  --  affine expression of JuMPeRVars and numbers.
"""
typealias AdaptExpr JuMP.GenericAffExpr{Float64,JuMPeRVar}
AdaptExpr() = zero(AdaptExpr)
Base.convert(::Type{AdaptExpr}, c::Number) =
    AdaptExpr(JuMPeRVar[ ], Float64[ ], 0.0)
Base.convert(::Type{AdaptExpr}, x::JuMPeRVar) =
    AdaptExpr(JuMPeRVar[x],Float64[1], 0.0)
Base.convert(::Type{AdaptExpr}, aff::AffExpr) =
    AdaptExpr(copy(aff.vars), copy(aff.coeffs), aff.constant)


"""
    AdaptConstraint

A constraint with just JuMPeRVars and numbers (i.e., `AdaptExpr`).
"""
typealias AdaptConstraint JuMP.GenericRangeConstraint{AdaptExpr}
function JuMP.addconstraint(m::Model, c::AdaptConstraint)
    rm = get_robust(m)::RobustModelExt
    push!(rm.adapt_constraints, c)
    return ConstraintRef{Model,AdaptConstraint}(m, length(rm.adapt_constraints))
end
