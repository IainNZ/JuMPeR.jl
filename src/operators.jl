#-----------------------------------------------------------------------
# JuMPeR  --  JuMP Extension for Robust Optimization
# http://github.com/IainNZ/JuMPeR.jl
#-----------------------------------------------------------------------
# Copyright (c) 2016: Iain Dunning
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
#-----------------------------------------------------------------------
# src/robustops.jl
# All the overloads for the new robust types introduced by JuMPeR.
# We share the same loose ordering as JuMP, just extended with the
# three new types:
# 1. Number
# 2. Variable
# 3. [Generic]Norm
# 4. [Generic]AffExpr   [Number * Variable]
# 5. QuadExpr <- We don't support any interactions with QuadExpr
# 6. [Generic]NormExpr
# 7. Uncertain
# 8. UncExpr            [Number * Uncertain]
# 9. UncVarExpr         [UncExpr * (Variable,Adaptive)]
# 10. Adaptive
# 11. AdaptExpr         [Number * Adaptive]
#-----------------------------------------------------------------------

# Define error message names for the sake of consistency
AFF = "an affine function of variables"
UNC = "an uncertain parameter"
UNE = "an affine function of uncertain parameters"
UVE = "an affine function of variables and uncertain parameters"
UNM = "a norm of uncertain parameters"
ADP = "an adaptive variable"


#-----------------------------------------------------------------------
# 1. Number
# Number--Uncertain
Base.:+(lhs::Number, rhs::Uncertain) = UncExpr([rhs], [  1], lhs)
Base.:-(lhs::Number, rhs::Uncertain) = UncExpr([rhs], [ -1], lhs)
Base.:*(lhs::Number, rhs::Uncertain) = UncExpr([rhs], [lhs], 0.0)
Base.:/(lhs::Number, rhs::Uncertain) = error("Cannot divide a number by $UNC")
# Number--UncExpr        - handled by JuMP
# Number--UncVarExpr     - handled by JuMP
# Number--Adaptive
Base.:+(c::Number, x::Adaptive) = AdaptExpr(Adaptive[x], Float64[+1], c)
Base.:-(c::Number, x::Adaptive) = AdaptExpr(Adaptive[x], Float64[-1], c)
Base.:*(c::Number, x::Adaptive) = AdaptExpr(Adaptive[x], Float64[ c], 0)
Base.:/(c::Number, x::Adaptive) = error("Cannot divide a number by $ADP")
# Number--AdaptExpr
Base.:+(c::Number, x::AdaptExpr) = AdaptExpr(copy(x.vars), copy(x.coeffs), c + x.constant)
Base.:-(c::Number, x::AdaptExpr) = AdaptExpr(copy(x.vars),     -x.coeffs , c - x.constant)
Base.:*(c::Number, x::AdaptExpr) = AdaptExpr(copy(x.vars),  c * x.coeffs , c * x.constant)
Base.:/(c::Number, x::AdaptExpr) = error("Cannot divide a number by $AFF")


#-----------------------------------------------------------------------
# 2. Variable
# Variable--Uncertain
Base.:+(lhs::Variable, rhs::Uncertain) = UncVarExpr([lhs],[UncExpr(  1)], UncExpr(rhs))
Base.:-(lhs::Variable, rhs::Uncertain) = UncVarExpr([lhs],[UncExpr(  1)],-UncExpr(rhs))
Base.:*(lhs::Variable, rhs::Uncertain) = UncVarExpr([lhs],[UncExpr(rhs)], UncExpr())
Base.:/(lhs::Variable, rhs::Uncertain) = error("Cannot divide a variable by $UNC")
# Variable--UncExpr
Base.:+(lhs::Variable, rhs::UncExpr) = UncVarExpr([lhs],[UncExpr(1)],       rhs)
Base.:-(lhs::Variable, rhs::UncExpr) = UncVarExpr([lhs],[UncExpr(1)],      -rhs)
Base.:*(lhs::Variable, rhs::UncExpr) = UncVarExpr([lhs],[       rhs], UncExpr())
Base.:/(lhs::Variable, rhs::UncExpr) = error("Cannot divide a variable by $UNE")
# Variable--UncVarExpr
Base.:+(lhs::Variable, rhs::UncVarExpr) = UncVarExpr(vcat(rhs.vars,lhs),vcat(      rhs.coeffs,UncExpr(1)), rhs.constant)
Base.:-(lhs::Variable, rhs::UncVarExpr) = UncVarExpr(vcat(rhs.vars,lhs),vcat(-1 .* rhs.coeffs,UncExpr(1)),-rhs.constant)
Base.:*(lhs::Variable, rhs::UncVarExpr) = error("Cannot multiply a variable by $UVE")
Base.:/(lhs::Variable, rhs::UncVarExpr) = error("Cannot divide a variable by $UVE")
# Variable--Adaptive
Base.:+(v::Variable, x::Adaptive) = AdaptExpr([v,x], [1,+1], 0)
Base.:-(v::Variable, x::Adaptive) = AdaptExpr([v,x], [1,-1], 0)
Base.:*(v::Variable, x::Adaptive) = error("Cannot multiply a variable by $ADP")
Base.:/(v::Variable, x::Adaptive) = error("Cannot divide a variable by $ADP")
# Variable--AdaptExpr
Base.:+(v::Variable, x::AdaptExpr) = AdaptExpr(vcat(v, x.vars), vcat(1,  x.coeffs),  x.constant)
Base.:-(v::Variable, x::AdaptExpr) = AdaptExpr(vcat(v, x.vars), vcat(1, -x.coeffs), -x.constant)
Base.:*(v::Variable, x::AdaptExpr) = error("Cannot multiply a variable by $AFF")
Base.:/(v::Variable, x::AdaptExpr) = error("Cannot divide a variable by $AFF")


#-----------------------------------------------------------------------
# 3. Norm
# Norm--Uncertain
Base.:*(lhs::GenericNorm{P,C,V}, rhs::Uncertain) where {P,C,V<:Uncertain} = error("Cannot multiply $UNM by $UNC")
# Norm--UncExpr
Base.:/(lhs::GenericNorm{P,C,V}, rhs::UncExpr) where {P,C,V<:Uncertain} = error("Cannot divide $UNM by $UNE")
# Norm--UncVarExpr
Base.:/(lhs::GenericNorm{P,C,V}, rhs::UncVarExpr) where {P,C,V<:Uncertain} = error("Cannot divide $UNM by $UVE")


#-----------------------------------------------------------------------
# 4. AffExpr
# AffExpr--Uncertain
Base.:+(lhs::AffExpr, rhs::Uncertain) = UncVarExpr(lhs.vars, map(UncExpr,lhs.coeffs), UncExpr([rhs],[ 1.],lhs.constant))
Base.:-(lhs::AffExpr, rhs::Uncertain) = UncVarExpr(lhs.vars, map(UncExpr,lhs.coeffs), UncExpr([rhs],[-1.],lhs.constant))
Base.:*(lhs::AffExpr, rhs::Uncertain) = UncVarExpr(lhs.vars,         rhs.*lhs.coeffs , UncExpr([rhs],[lhs.constant],0.0))
Base.:/(lhs::AffExpr, rhs::Uncertain) = error("Cannot divide $AFF by $UNC")
# AffExpr--UncExpr
Base.:+(lhs::AffExpr, rhs::UncExpr)  = UncVarExpr(lhs.vars, map(UncExpr,lhs.coeffs), lhs.constant+rhs)
Base.:-(lhs::AffExpr, rhs::UncExpr)  = UncVarExpr(lhs.vars, map(UncExpr,lhs.coeffs), lhs.constant-rhs)
Base.:*(lhs::AffExpr, rhs::UncExpr)  = UncVarExpr(lhs.vars,         rhs.*lhs.coeffs , lhs.constant*rhs)
Base.:/(lhs::AffExpr, rhs::UncExpr)  = error("Cannot divide $AFF by $UNE")
# AffExpr--UncVarExpr
Base.:+(lhs::AffExpr, rhs::UncVarExpr) = UncVarExpr(
  vcat(lhs.vars, rhs.vars),
  vcat(map(UncExpr,lhs.coeffs), rhs.coeffs),
  lhs.constant + rhs.constant)
Base.:-(lhs::AffExpr, rhs::UncVarExpr) = UncVarExpr(
  vcat(lhs.vars, rhs.vars),
  vcat(map(UncExpr,lhs.coeffs), -1 .* rhs.coeffs),
  lhs.constant - rhs.constant)
Base.:*(lhs::AffExpr, rhs::UncVarExpr) = error("Cannot multiply $AFF by $UVE")
Base.:/(lhs::AffExpr, rhs::UncVarExpr) = error("Cannot divide $AFF by $UVE")
# AffExpr--GenericNormExpr{Uncertain}
Base.:+(lhs::AffExpr, rhs::GenericNormExpr{P,Float64,Uncertain}) where P =
    length(lhs.vars) == 0 ? lhs.constant + rhs : error("Cannot add $AFF by $UNM")
Base.:-(lhs::AffExpr, rhs::GenericNormExpr{P,Float64,Uncertain}) where P =
    length(lhs.vars) == 0 ? lhs.constant - rhs : error("Cannot substract $AFF by $UNM")
# AffExpr--Adaptive
Base.:+(a::AffExpr, x::Adaptive) = AdaptExpr(vcat(a.vars, x), vcat(a.coeffs,  1), a.constant)
Base.:-(a::AffExpr, x::Adaptive) = AdaptExpr(vcat(a.vars, x), vcat(a.coeffs, -1), a.constant)
Base.:*(a::AffExpr, x::Adaptive) = error("Cannot multiply $AFF by $ADP")
Base.:/(a::AffExpr, x::Adaptive) = error("Cannot divide $AFF by $ADP")
# AffExpr--AdaptiveAffExpr
Base.:+(a::AffExpr, b::AdaptExpr) = AdaptExpr(vcat(a.vars,   b.vars),
                                                vcat(a.coeffs, b.coeffs),
                                                a.constant + b.constant)
Base.:-(a::AffExpr, b::AdaptExpr) = AdaptExpr(vcat(a.vars,   b.vars),
                                                vcat(a.coeffs,-b.coeffs),
                                                a.constant - b.constant)
Base.:*(a::AffExpr, b::AdaptExpr) = error("Cannot multiply $AFF by $AFF")
Base.:/(a::AffExpr, b::AdaptExpr) = error("Cannot divide $AFF by $AFF")


#-----------------------------------------------------------------------
# 5. QuadExpr
# Nothing additional supported


#-----------------------------------------------------------------------
# 6. GenericNormExpr
# GenericNormExpr--Uncertain
Base.:+(lhs::GenericNormExpr{P,C,V},rhs::Uncertain) where {P, C, V <: Uncertain} = error("Cannot add $UNM to $UNC")
Base.:-(lhs::GenericNormExpr{P,C,V},rhs::Uncertain) where {P, C, V <: Uncertain} = error("Cannot subtract $UNM by $UNC")
Base.:*(lhs::GenericNormExpr{P,C,V},rhs::Uncertain) where {P, C, V <: Uncertain} = error("Cannot multiply $UNM by $UNC")
Base.:/(lhs::GenericNormExpr{P,C,V},rhs::Uncertain) where {P, C, V <: Uncertain} = error("Cannot divide $UNM by $UNC")
# GenericNormExpr--UncExpr
Base.:+(lhs::GenericNormExpr{P,C,V},rhs::UncExpr) where {P, C, V <: Uncertain} = error("Cannot add $UNM to $UNE")
Base.:-(lhs::GenericNormExpr{P,C,V},rhs::UncExpr) where {P, C, V <: Uncertain} = error("Cannot subtract $UNM by $UNE")
Base.:*(lhs::GenericNormExpr{P,C,V},rhs::UncExpr) where {P, C, V <: Uncertain} = error("Cannot multiply $UNM by $UNE")
Base.:/(lhs::GenericNormExpr{P,C,V},rhs::UncExpr) where {P, C, V <: Uncertain} = error("Cannot divide $UNM by $UNE")
# GenericNormExpr--UncVarExpr
Base.:+(lhs::GenericNormExpr{P,C,V},rhs::UncVarExpr) where {P, C, V <: Uncertain} = error("Cannot add $UNM to $UVE")
Base.:-(lhs::GenericNormExpr{P,C,V},rhs::UncVarExpr) where {P, C, V <: Uncertain} = error("Cannot subtract $UNM by $UVE")
Base.:*(lhs::GenericNormExpr{P,C,V},rhs::UncVarExpr) where {P, C, V <: Uncertain} = error("Cannot multiply $UNM by $UVE")
Base.:/(lhs::GenericNormExpr{P,C,V},rhs::UncVarExpr) where {P, C, V <: Uncertain} = error("Cannot divide $UNM by $UVE")


#-----------------------------------------------------------------------
# 7. Uncertain
Base.:-(lhs::Uncertain) = UncExpr(-1,lhs)
# Uncertain--Number
Base.:+(lhs::Uncertain, rhs::Number) = UncExpr([lhs],[1.0], rhs)
Base.:-(lhs::Uncertain, rhs::Number) = UncExpr([lhs],[1.0],-rhs)
Base.:*(lhs::Uncertain, rhs::Number) = UncExpr(  rhs, lhs)
Base.:/(lhs::Uncertain, rhs::Number) = UncExpr(1/rhs, lhs)
# Uncertain--Variable
Base.:+(lhs::Uncertain, rhs::Variable) = Base.:+(rhs, lhs)
Base.:-(lhs::Uncertain, rhs::Variable) = UncVarExpr([rhs],[UncExpr(-1)],UncExpr(lhs))
Base.:*(lhs::Uncertain, rhs::Variable) = Base.:*(rhs, lhs)
Base.:/(lhs::Uncertain, rhs::Variable) = error("Cannot divide $UNC by a variable")
# Uncertain--GenericNorm
Base.:*(lhs::Uncertain, rhs::GenericNorm{P,C,V}) where {P, C, V <: Uncertain} = error("Cannot multiply $UNC by $UNM")
Base.:/(lhs::Uncertain, rhs::GenericNorm{P,C,V}) where {P, C, V <: Uncertain} = error("Cannot multiply $UNC by $UNM")
# Uncertain--AffExpr
Base.:+(lhs::Uncertain, rhs::AffExpr) = Base.:+(rhs, lhs)
Base.:-(lhs::Uncertain, rhs::AffExpr) = UncVarExpr(rhs.vars, map(UncExpr,-rhs.coeffs), UncExpr([lhs],[1.0],-rhs.constant))
Base.:*(lhs::Uncertain, rhs::AffExpr) = Base.:*(rhs, lhs)
Base.:/(lhs::Uncertain, rhs::AffExpr) = error("Cannot divide $UNC by $AFF")
# Uncertain--GenericNormExpr
Base.:+(lhs::Uncertain, rhs::GenericNormExpr{P,C,V}) where {P, C, V <: Uncertain} = error("Cannot add $UNC to $UNM")
Base.:-(lhs::Uncertain, rhs::GenericNormExpr{P,C,V}) where {P, C, V <: Uncertain} = error("Cannot subtract $UNC by $UNM")
Base.:*(lhs::Uncertain, rhs::GenericNormExpr{P,C,V}) where {P, C, V <: Uncertain} = error("Cannot multiply $UNC by $UNM")
Base.:/(lhs::Uncertain, rhs::GenericNormExpr{P,C,V}) where {P, C, V <: Uncertain} = error("Cannot divide $UNC by $UNM")
# Uncertain--Uncertain
Base.:+(lhs::Uncertain, rhs::Uncertain) = UncExpr([lhs,rhs], Float64[1, 1], 0.0)
Base.:-(lhs::Uncertain, rhs::Uncertain) = UncExpr([lhs,rhs], Float64[1,-1], 0.0)
Base.:*(lhs::Uncertain, rhs::Uncertain) = error("Cannot multiply $UNC by $UNC")
Base.:/(lhs::Uncertain, rhs::Uncertain) = error("Cannot divide $UNC by $UNC")
# Uncertain--UncExpr
Base.:+(lhs::Uncertain, rhs::UncExpr) = UncExpr(vcat(lhs,rhs.vars),vcat(1.0, rhs.coeffs), rhs.constant)
Base.:-(lhs::Uncertain, rhs::UncExpr) = UncExpr(vcat(lhs,rhs.vars),vcat(1.0,-rhs.coeffs),-rhs.constant)
Base.:*(lhs::Uncertain, rhs::UncExpr) = error("Cannot multiply $UNC by $UNE")
Base.:/(lhs::Uncertain, rhs::UncExpr) = error("Cannot divide $UNC by $UNE")
# Uncertain--UncVarExpr
Base.:+(lhs::Uncertain, rhs::UncVarExpr) = UncVarExpr(rhs.vars,    rhs.coeffs,lhs+rhs.constant)
Base.:-(lhs::Uncertain, rhs::UncVarExpr) = UncVarExpr(rhs.vars,-1 .* rhs.coeffs,lhs-rhs.constant)
Base.:*(lhs::Uncertain, rhs::UncVarExpr) = error("Cannot multiply $UNC by $UVE")
Base.:/(lhs::Uncertain, rhs::UncVarExpr) = error("Cannot divide $UNC by $UVE")
# Uncertain--Adaptive
Base.:+(u::Uncertain, x::Adaptive) = UncVarExpr([x], [ 1], u)
Base.:-(u::Uncertain, x::Adaptive) = UncVarExpr([x], [-1], u)
Base.:*(u::Uncertain, x::Adaptive) = UncVarExpr([x], [ u], 0)
Base.:/(u::Uncertain, x::Adaptive) = error("Cannot divide $UNC by $ADP")
# Uncertain--AdaptExpr
Base.:+(u::Uncertain, x::AdaptExpr) = UncVarExpr(copy(x.vars), copy(x.coeffs), u + x.constant)
Base.:-(u::Uncertain, x::AdaptExpr) = UncVarExpr(copy(x.vars),     -x.coeffs , u - x.constant)
Base.:*(u::Uncertain, x::AdaptExpr) = UncVarExpr(copy(x.vars), u .* x.coeffs , u * x.constant)
Base.:/(u::Uncertain, x::AdaptExpr) = error("Cannout divide $UNC by $AFF")


#-----------------------------------------------------------------------
# 8. UncExpr
# UncExpr--Number        - handled by JuMP
# UncExpr--Variable
Base.:+(lhs::UncExpr, rhs::Variable) = Base.:+(rhs,lhs)
Base.:-(lhs::UncExpr, rhs::Variable) = UncVarExpr([rhs],[UncExpr(-1)],lhs)
Base.:*(lhs::UncExpr, rhs::Variable) = Base.:*(rhs,lhs)
Base.:/(lhs::UncExpr, rhs::Variable) = error("Cannot divide $UNE by a variable")
# UncExpr--AffExpr
Base.:+(lhs::UncExpr, rhs::AffExpr) = Base.:+( rhs,lhs)
Base.:-(lhs::UncExpr, rhs::AffExpr) = Base.:+(-rhs,lhs)
Base.:*(lhs::UncExpr, rhs::AffExpr) = Base.:*( rhs,lhs)
Base.:/(lhs::UncExpr, rhs::AffExpr) = error("Cannot divide $UNE by $AFF")
# UncExpr--GenericNormExpr
Base.:+(lhs::UncExpr, rhs::GenericNormExpr{P,C,V}) where {P, C, V <: Uncertain} = error("Cannot add $UNE to $UNM")
Base.:-(lhs::UncExpr, rhs::GenericNormExpr{P,C,V}) where {P, C, V <: Uncertain} = error("Cannot subtract $UNE by $UNM")
Base.:*(lhs::UncExpr, rhs::GenericNormExpr{P,C,V}) where {P, C, V <: Uncertain} = error("Cannot multiply $UNE by $UNM")
Base.:/(lhs::UncExpr, rhs::GenericNormExpr{P,C,V}) where {P, C, V <: Uncertain} = error("Cannot divide $UNE by $UNM")
# UncExpr--Uncertain
Base.:+(lhs::UncExpr, rhs::Uncertain) = Base.:+(rhs,lhs)
Base.:-(lhs::UncExpr, rhs::Uncertain) = UncExpr(vcat(rhs,lhs.vars),vcat(-1.0,lhs.coeffs),lhs.constant)
Base.:*(lhs::UncExpr, rhs::Uncertain) = Base.:*(rhs,lhs)
Base.:/(lhs::UncExpr, rhs::Uncertain) = error("Cannot divide $UNE by $UNC")
# UncExpr--UncExpr
Base.:*(lhs::UncExpr, rhs::UncExpr) = error("Cannot multiply $UNE by $UNE")
Base.:/(lhs::UncExpr, rhs::UncExpr) = error("Cannot divide $UNE by $UNE")
# UncExpr--UncVarExpr
Base.:+(lhs::UncExpr, rhs::UncVarExpr) = UncVarExpr(rhs.vars,      rhs.coeffs,lhs+rhs.constant)
Base.:-(lhs::UncExpr, rhs::UncVarExpr) = UncVarExpr(rhs.vars,-1 .* rhs.coeffs,lhs-rhs.constant)
Base.:*(lhs::UncExpr, rhs::UncVarExpr) = (length(lhs.vars) == 0) ?
                                        lhs.constant * rhs : # LHS is just a constant, so OK
                                        error("Cannot multiply $UNE by $UVE")
Base.:/(lhs::UncExpr, rhs::UncVarExpr) = error("Cannot divide $UNE by $UVE")
# UncExpr--Adaptive
Base.:+(u::UncExpr, x::Adaptive) = UncVarExpr([x], [ 1], u)
Base.:-(u::UncExpr, x::Adaptive) = UncVarExpr([x], [-1], u)
Base.:*(u::UncExpr, x::Adaptive) = UncVarExpr([x], [ u], 0)
Base.:/(u::UncExpr, x::Adaptive) = error("Cannot divide $UNE by $ADP")
# UncExpr--AdaptExpr
Base.:+(u::UncExpr, x::AdaptExpr) = UncVarExpr(copy(x.vars), copy(x.coeffs), u + x.constant)
Base.:-(u::UncExpr, x::AdaptExpr) = UncVarExpr(copy(x.vars),     -x.coeffs , u - x.constant)
Base.:*(u::UncExpr, x::AdaptExpr) = UncVarExpr(copy(x.vars), u .* x.coeffs , u * x.constant)
Base.:/(u::UncExpr, x::AdaptExpr) = error("Cannot divide $UNE by $AFF")


#-----------------------------------------------------------------------
# 9. UncVarExpr
# UncVarExpr--Number     - handled by JuMP
# UncVarExpr--Variable
Base.:+(lhs::UncVarExpr, rhs::Variable) = Base.:+(rhs,lhs)
Base.:-(lhs::UncVarExpr, rhs::Variable) = UncVarExpr(vcat(lhs.vars,rhs),vcat(lhs.coeffs,UncExpr(-1)), lhs.constant)
Base.:*(lhs::UncVarExpr, rhs::Variable) = error("Cannot multiply $UVE by a variable")
Base.:/(lhs::UncVarExpr, rhs::Variable) = error("Cannot divide $UVE by a variable")
# UncVarExpr--AffExpr
Base.:+(lhs::UncVarExpr, rhs::AffExpr) = Base.:+(rhs,lhs)
Base.:-(lhs::UncVarExpr, rhs::AffExpr) = UncVarExpr(
  vcat(lhs.vars,    rhs.vars),
  vcat(lhs.coeffs, -rhs.coeffs),
  lhs.constant - rhs.constant)
Base.:*(lhs::UncVarExpr, rhs::AffExpr) = error("Cannot multiply $UVE by $AFF")
Base.:/(lhs::UncVarExpr, rhs::AffExpr) = error("Cannot divide $UVE by $AFF")
# UncVarExpr--GenericNormExpr
Base.:+(lhs::UncVarExpr, rhs::GenericNormExpr{P,C,V}) where {P, C, V <: Uncertain} = error("Cannot add $UVE to $UNM")
Base.:-(lhs::UncVarExpr, rhs::GenericNormExpr{P,C,V}) where {P, C, V <: Uncertain} = error("Cannot subtract $UVE by $UNM")
Base.:*(lhs::UncVarExpr, rhs::GenericNormExpr{P,C,V}) where {P, C, V <: Uncertain} = error("Cannot multiply $UVE by $UNM")
Base.:/(lhs::UncVarExpr, rhs::GenericNormExpr{P,C,V}) where {P, C, V <: Uncertain} = error("Cannot divide $UVE by $UNM")
# UncVarExpr--Uncertain
Base.:+(lhs::UncVarExpr, rhs::Uncertain) = UncVarExpr(lhs.vars,lhs.coeffs,lhs.constant+rhs)
Base.:-(lhs::UncVarExpr, rhs::Uncertain) = UncVarExpr(lhs.vars,lhs.coeffs,lhs.constant-rhs)
Base.:*(lhs::UncVarExpr, rhs::Uncertain) = error("Cannot multiply $UVE by $UNC")
Base.:/(lhs::UncVarExpr, rhs::Uncertain) = error("Cannot divide $UVE by $UNC")
# UncVarExpr--UncExpr
Base.:+(lhs::UncVarExpr, rhs::UncExpr) = UncVarExpr(lhs.vars,lhs.coeffs,lhs.constant+rhs)
Base.:-(lhs::UncVarExpr, rhs::UncExpr) = UncVarExpr(lhs.vars,lhs.coeffs,lhs.constant-rhs)
Base.:*(lhs::UncVarExpr, rhs::UncExpr) = error("Cannot multiply $UVE by $UNE")
Base.:/(lhs::UncVarExpr, rhs::UncExpr) = error("Cannot divide $UVE by $UNE")
# UncVarExpr--UncVarExpr
Base.:*(lhs::UncVarExpr, rhs::UncVarExpr) = error("Cannot multiply $UVE by $UVE")
Base.:/(lhs::UncVarExpr, rhs::UncVarExpr) = error("Cannot divide $UVE by $UVE")
# UncVarExpr--Adaptive
Base.:+(a::UncVarExpr, x::Adaptive) = UncVarExpr(vcat(a.vars, x),
                                           vcat(a.coeffs, UncExpr(1)),
                                           a.constant)
Base.:-(a::UncVarExpr, x::Adaptive) = UncVarExpr(vcat(a.vars, x),
                                           vcat(a.coeffs, UncExpr(-1)),
                                           a.constant)
Base.:*(a::UncVarExpr, x::Adaptive) = error("Cannot multiply $UVE by $ADP")
Base.:/(a::UncVarExpr, x::Adaptive) = error("Cannot divide $UVE by $ADP")
# UncVarExpr--AdaptExpr
Base.:+(a::UncVarExpr, b::AdaptExpr) = UncVarExpr(vcat(a.vars, b.vars),
                                               vcat(a.coeffs, map(UncExpr, b.coeffs)),
                                               a.constant + b.constant)
Base.:-(a::UncVarExpr, b::AdaptExpr) = UncVarExpr(vcat(a.vars, b.vars),
                                               vcat(a.coeffs, map(UncExpr,-b.coeffs)),
                                               a.constant - b.constant)
Base.:*(a::UncVarExpr, b::AdaptExpr) = error("Cannot multiply $UVE by $AFF")
Base.:/(a::UncVarExpr, b::AdaptExpr) = error("Cannot divide $UVE by $AFF")


#-----------------------------------------------------------------------
# 10. Adaptive
# Adaptive--Number
Base.:+(x::Adaptive, c::Number) = +(  c, x)
Base.:-(x::Adaptive, c::Number) = +( -c, x)
Base.:*(x::Adaptive, c::Number) = *(  c, x)
Base.:/(x::Adaptive, c::Number) = *(1/c, x)
# Adaptive--Variable
Base.:+(x::Adaptive, v::Variable) = AdaptExpr([x,v], [1,+1], 0)
Base.:-(x::Adaptive, v::Variable) = AdaptExpr([x,v], [1,-1], 0)
Base.:*(x::Adaptive, v::Variable) = error("Cannot multiply $ADP by a variable")
Base.:/(x::Adaptive, v::Variable) = error("Cannot divide $ADP by a variable")
# Adaptive--AffExpr
Base.:+(x::Adaptive, a::AffExpr) = AdaptExpr(vcat(x, a.vars), vcat(1,  a.coeffs),  a.constant)
Base.:-(x::Adaptive, a::AffExpr) = AdaptExpr(vcat(x, a.vars), vcat(1, -a.coeffs), -a.constant)
Base.:*(x::Adaptive, a::AffExpr) = error("Cannot multiply $ADP by $AFF")
Base.:/(x::Adaptive, a::AffExpr) = error("Cannot divide $ADP by $AFF")
# Adaptive--Uncertain
Base.:+(x::Adaptive, u::Uncertain) = UncVarExpr([x], [1],  u)
Base.:-(x::Adaptive, u::Uncertain) = UncVarExpr([x], [1], -u)
Base.:*(x::Adaptive, u::Uncertain) = UncVarExpr([x], [u],  0)
Base.:/(x::Adaptive, u::Uncertain) = error("Cannot divide $ADP by $UNC")
# Adaptive--UncExpr
Base.:+(x::Adaptive, u::UncExpr) = UncVarExpr([x], [1],  u)
Base.:-(x::Adaptive, u::UncExpr) = UncVarExpr([x], [1], -u)
Base.:*(x::Adaptive, u::UncExpr) = UncVarExpr([x], [u],  0)
Base.:/(x::Adaptive, u::UncExpr) = error("Cannot divide $ADP by $UNE")
# Adaptive--UncVarExpr
Base.:+(x::Adaptive, u::UncVarExpr) = UncVarExpr(vcat(x, u.vars), vcat(1,  u.coeffs),  u.constant)
Base.:-(x::Adaptive, u::UncVarExpr) = UncVarExpr(vcat(x, u.vars), vcat(1, -u.coeffs), -u.constant)
Base.:*(x::Adaptive, u::UncVarExpr) = error("Cannot multiply $ADP by $UVE")
Base.:/(x::Adaptive, u::UncVarExpr) = error("Cannot divide $ADP by $UVE")
# Adaptive--Adaptive
Base.:+(a::Adaptive, b::Adaptive) = AdaptExpr(Adaptive[a,b], Float64[1, 1], 0)
Base.:-(a::Adaptive, b::Adaptive) = AdaptExpr(Adaptive[a,b], Float64[1,-1], 0)
Base.:*(a::Adaptive, b::Adaptive) = error("Cannot multiply $ADP by $ADP")
Base.:/(a::Adaptive, b::Adaptive) = error("Cannot divide $ADP by $ADP")
# Adaptive--AdaptExpr
Base.:+(a::Adaptive, b::AdaptExpr) = AdaptExpr(vcat(a, b.vars), vcat(1,  b.coeffs),  b.constant)
Base.:-(a::Adaptive, b::AdaptExpr) = AdaptExpr(vcat(a, b.vars), vcat(1, -b.coeffs), -b.constant)
Base.:*(a::Adaptive, b::AdaptExpr) = error("Cannot multiply $ADP by $AFF")
Base.:/(a::Adaptive, b::AdaptExpr) = error("Cannot divide $ADP by $AFF")


#-----------------------------------------------------------------------
# 11. AdaptExpr
# AdaptExpr--Number
Base.:+(a::AdaptExpr, c::Number) = AdaptExpr(copy(a.vars), copy(a.coeffs), a.constant + c)
Base.:-(a::AdaptExpr, c::Number) = AdaptExpr(copy(a.vars), copy(a.coeffs), a.constant - c)
Base.:*(a::AdaptExpr, c::Number) = AdaptExpr(copy(a.vars), copy(a.coeffs)*c, a.constant*c)
Base.:/(a::AdaptExpr, c::Number) = AdaptExpr(copy(a.vars), copy(a.coeffs)/c, a.constant/c)
# AdaptExpr--Variable
Base.:+(x::AdaptExpr, v::Variable) = AdaptExpr(vcat(x.vars, v), vcat(x.coeffs,  1), x.constant)
Base.:-(x::AdaptExpr, v::Variable) = AdaptExpr(vcat(x.vars, v), vcat(x.coeffs, -1), x.constant)
Base.:*(x::AdaptExpr, v::Variable) = error("Cannot multiply $AFF by a variable")
Base.:/(x::AdaptExpr, v::Variable) = error("Cannot divide $AFF by a variable")
# AdaptExpr--AffExpr
Base.:+(x::AdaptExpr, a::AffExpr) = AdaptExpr(vcat(x.vars, a.vars), vcat(x.coeffs,  a.coeffs), x.constant + a.constant)
Base.:-(x::AdaptExpr, a::AffExpr) = AdaptExpr(vcat(x.vars, a.vars), vcat(x.coeffs, -a.coeffs), x.constant - a.constant)
Base.:*(x::AdaptExpr, a::AffExpr) = error("Cannot multiply $AFF by $AFF")
Base.:/(x::AdaptExpr, a::AffExpr) = error("Cannot divide $AFF by $AFF")
# AdaptExpr--Uncertain
Base.:+(a::AdaptExpr, u::Uncertain) = UncVarExpr(copy(a.vars), copy(a.coeffs), a.constant + u)
Base.:-(a::AdaptExpr, u::Uncertain) = UncVarExpr(copy(a.vars), copy(a.coeffs), a.constant - u)
Base.:*(a::AdaptExpr, u::Uncertain) = UncVarExpr(copy(a.vars), copy(a.coeffs) .* u, a.constant * u)
Base.:/(a::AdaptExpr, u::Uncertain) = error("Cannot divide $AFF by $UNC")
# AdaptExpr--UncExpr
Base.:+(a::AdaptExpr, b::UncExpr) = +( b, a)
Base.:-(a::AdaptExpr, b::UncExpr) = +(-b, a)
Base.:*(a::AdaptExpr, b::UncExpr) = *( b, a)
Base.:/(a::AdaptExpr, b::UncExpr) = /( b, a)  # not quite
# AdaptExpr--UncVarExpr
Base.:+(a::AdaptExpr, b::UncVarExpr) = UncVarExpr(vcat(a.vars, b.vars),
                                                 vcat(map(UncExpr,a.coeffs), b.coeffs),
                                                 a.constant + b.constant)
Base.:-(a::AdaptExpr, b::UncVarExpr) = UncVarExpr(vcat(a.vars, b.vars),
                                                 vcat(map(UncExpr,a.coeffs), -b.coeffs),
                                                 a.constant - b.constant)
Base.:*(a::AdaptExpr, b::UncVarExpr) = error("Cannot multiply $AFF by $UVE")
Base.:/(a::AdaptExpr, b::UncVarExpr) = error("Cannot divide $AFF by $UVE")
# AdaptExpr--Adaptive
Base.:+(a::AdaptExpr, b::Adaptive) = AdaptExpr(vcat(a.vars,b), vcat(a.coeffs, 1), a.constant)
Base.:-(a::AdaptExpr, b::Adaptive) = AdaptExpr(vcat(a.vars,b), vcat(a.coeffs,-1), a.constant)
Base.:*(a::AdaptExpr, b::Adaptive) = error("Cannot multiply $AFF by $ADP")
Base.:/(a::AdaptExpr, b::Adaptive) = error("Cannot divide $AFF by $ADP")
# AdaptExpr--AdaptExpr
Base.:+(a::AdaptExpr, b::AdaptExpr) = AdaptExpr(vcat(a, b.vars), vcat(1,  b.coeffs), a.constant+b.constant)
Base.:-(a::AdaptExpr, b::AdaptExpr) = AdaptExpr(vcat(a, b.vars), vcat(1, -b.coeffs), a.constant-b.constant)
Base.:*(a::AdaptExpr, b::AdaptExpr) = error("Cannot multiply $AFF by $AFF")
Base.:/(a::AdaptExpr, b::AdaptExpr) = error("Cannot divide $AFF by $AFF")


#-----------------------------------------------------------------------
# For matrix operations
Base.promote_rule(::Type{R}, ::Type{Uncertain})  where R<:Real  = UncExpr
Base.promote_rule(::Type{R}, ::Type{UncExpr})    where R<:Real  = UncExpr
Base.promote_rule(::Type{R}, ::Type{UncVarExpr}) where R<:Real  = UncVarExpr
Base.promote_rule(::Type{R}, ::Type{Adaptive})   where R<:Real  = AdaptExpr
Base.promote_rule(::Type{R}, ::Type{AdaptExpr})  where R<:Real  = AdaptExpr

Base.promote_rule(::Type{Uncertain}, ::Type{R}) where R<:Real   = UncExpr
Base.promote_rule(::Type{Uncertain}, ::Type{Uncertain})         = UncExpr
Base.promote_rule(::Type{Uncertain}, ::Type{UncExpr})           = UncExpr
Base.promote_rule(::Type{Uncertain}, ::Type{UncVarExpr})        = UncVarExpr
Base.promote_rule(::Type{Uncertain}, ::Type{Adaptive})          = UncVarExpr
Base.promote_rule(::Type{Uncertain}, ::Type{AdaptExpr})         = UncVarExpr

Base.promote_rule(::Type{UncExpr}, ::Type{R}) where R<:Real     = UncExpr
Base.promote_rule(::Type{UncExpr}, ::Type{Uncertain})           = UncExpr
Base.promote_rule(::Type{UncExpr}, ::Type{UncExpr})             = UncExpr
Base.promote_rule(::Type{UncExpr}, ::Type{UncVarExpr})          = UncVarExpr
Base.promote_rule(::Type{UncExpr}, ::Type{Adaptive})            = UncVarExpr
Base.promote_rule(::Type{UncExpr}, ::Type{AdaptExpr})           = UncVarExpr

Base.promote_rule(::Type{UncVarExpr}, ::Type{R}) where R<:Real  = UncVarExpr
Base.promote_rule(::Type{UncVarExpr}, ::Type{Uncertain})        = UncVarExpr
Base.promote_rule(::Type{UncVarExpr}, ::Type{UncExpr})          = UncVarExpr
Base.promote_rule(::Type{UncVarExpr}, ::Type{UncVarExpr})       = UncVarExpr
Base.promote_rule(::Type{UncVarExpr}, ::Type{Adaptive})         = UncVarExpr
Base.promote_rule(::Type{UncVarExpr}, ::Type{AdaptExpr})        = UncVarExpr

Base.promote_rule(::Type{Adaptive}, ::Type{R}) where R<:Real    = AdaptExpr
Base.promote_rule(::Type{Adaptive}, ::Type{Uncertain})          = UncVarExpr
Base.promote_rule(::Type{Adaptive}, ::Type{UncExpr})            = UncVarExpr
Base.promote_rule(::Type{Adaptive}, ::Type{UncVarExpr})         = UncVarExpr
Base.promote_rule(::Type{Adaptive}, ::Type{Adaptive})           = AdaptExpr
Base.promote_rule(::Type{Adaptive}, ::Type{AdaptExpr})          = AdaptExpr

Base.promote_rule(::Type{AdaptExpr}, ::Type{R}) where R<:Real   = AdaptExpr
Base.promote_rule(::Type{AdaptExpr}, ::Type{Uncertain})         = UncVarExpr
Base.promote_rule(::Type{AdaptExpr}, ::Type{UncExpr})           = UncVarExpr
Base.promote_rule(::Type{AdaptExpr}, ::Type{UncVarExpr})        = UncVarExpr
Base.promote_rule(::Type{AdaptExpr}, ::Type{Adaptive})          = AdaptExpr
Base.promote_rule(::Type{AdaptExpr}, ::Type{AdaptExpr})         = AdaptExpr
