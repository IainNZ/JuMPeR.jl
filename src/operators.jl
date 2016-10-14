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
const AFF = "an affine function of variables"
const UNC = "an uncertain parameter"
const UNE = "an affine function of uncertain parameters"
const UVE = "an affine function of variables and uncertain parameters"
const UNM = "a norm of uncertain parameters"
const ADP = "an adaptive variable"


#-----------------------------------------------------------------------
# 1. Number
# Number--Uncertain
(+)(lhs::Number, rhs::Uncertain) = UncExpr([rhs], [  1], lhs)
(-)(lhs::Number, rhs::Uncertain) = UncExpr([rhs], [ -1], lhs)
(*)(lhs::Number, rhs::Uncertain) = UncExpr(lhs, rhs)
(/)(lhs::Number, rhs::Uncertain) = error("Cannot divide a number by $UNC")
# Number--UncExpr        - handled by JuMP
# Number--UncVarExpr     - handled by JuMP
# Number--Adaptive
(+)(c::Number, x::Adaptive) = AdaptExpr(Adaptive[x], Float64[+1], c)
(-)(c::Number, x::Adaptive) = AdaptExpr(Adaptive[x], Float64[-1], c)
(*)(c::Number, x::Adaptive) = AdaptExpr(Adaptive[x], Float64[ c], 0)
(/)(c::Number, x::Adaptive) = error("Cannot divide a number by $ADP")
# Number--AdaptExpr
(+)(c::Number, x::AdaptExpr) = AdaptExpr(copy(x.vars), copy(x.coeffs), c + x.constant)
(-)(c::Number, x::AdaptExpr) = AdaptExpr(copy(x.vars),     -x.coeffs , c - x.constant)
(*)(c::Number, x::AdaptExpr) = AdaptExpr(copy(x.vars),  c * x.coeffs , c * x.constant)
(/)(c::Number, x::AdaptExpr) = error("Cannot divide a number by $AFF")


#-----------------------------------------------------------------------
# 2. Variable
# Variable--Uncertain
(+)(lhs::Variable, rhs::Uncertain) = UncVarExpr([lhs],[UncExpr(  1)], UncExpr(rhs))
(-)(lhs::Variable, rhs::Uncertain) = UncVarExpr([lhs],[UncExpr(  1)],-UncExpr(rhs))
(*)(lhs::Variable, rhs::Uncertain) = UncVarExpr([lhs],[UncExpr(rhs)], UncExpr())
(/)(lhs::Variable, rhs::Uncertain) = error("Cannot divide a variable by $UNC")
# Variable--UncExpr
(+)(lhs::Variable, rhs::UncExpr) = UncVarExpr([lhs],[UncExpr(1)],       rhs)
(-)(lhs::Variable, rhs::UncExpr) = UncVarExpr([lhs],[UncExpr(1)],      -rhs)
(*)(lhs::Variable, rhs::UncExpr) = UncVarExpr([lhs],[        rhs],UncExpr())
(/)(lhs::Variable, rhs::UncExpr) = error("Cannot divide a variable by $UNE")
# Variable--UncVarExpr
(+)(lhs::Variable, rhs::UncVarExpr) = UncVarExpr(vcat(rhs.vars,lhs),vcat(    rhs.coeffs,UncExpr(1)), rhs.constant)
(-)(lhs::Variable, rhs::UncVarExpr) = UncVarExpr(vcat(rhs.vars,lhs),vcat(-1.*rhs.coeffs,UncExpr(1)),-rhs.constant)
(*)(lhs::Variable, rhs::UncVarExpr) = error("Cannot multiply a variable by $UVE")
(/)(lhs::Variable, rhs::UncVarExpr) = error("Cannot divide a variable by $UVE")
# Variable--Adaptive
(+)(v::Variable, x::Adaptive) = AdaptExpr([v,x], [1,+1], 0)
(-)(v::Variable, x::Adaptive) = AdaptExpr([v,x], [1,-1], 0)
(*)(v::Variable, x::Adaptive) = error("Cannot multiply a variable by $ADP")
(/)(v::Variable, x::Adaptive) = error("Cannot divide a variable by $ADP")
# Variable--AdaptExpr
(+)(v::Variable, x::AdaptExpr) = AdaptExpr(vcat(v, x.vars), vcat(1,  x.coeffs),  x.constant)
(-)(v::Variable, x::AdaptExpr) = AdaptExpr(vcat(v, x.vars), vcat(1, -x.coeffs), -x.constant)
(*)(v::Variable, x::AdaptExpr) = error("Cannot multiply a variable by $AFF")
(/)(v::Variable, x::AdaptExpr) = error("Cannot divide a variable by $AFF")


#-----------------------------------------------------------------------
# 3. Norm
# Norm--Uncertain
(*){P,C,V<:Uncertain}(lhs::GenericNorm{P,C,V}, rhs::Uncertain) = error("Cannot multiply $UNM by $UNC")
# Norm--UncExpr
(/){P,C,V<:Uncertain}(lhs::GenericNorm{P,C,V}, rhs::UncExpr) = error("Cannot divide $UNM by $UNE")
# Norm--UncVarExpr
(/){P,C,V<:Uncertain}(lhs::GenericNorm{P,C,V}, rhs::UncVarExpr) = error("Cannot divide $UNM by $UVE")


#-----------------------------------------------------------------------
# 4. AffExpr
# AffExpr--Uncertain
(+)(lhs::AffExpr, rhs::Uncertain) = UncVarExpr(lhs.vars, map(UncExpr,lhs.coeffs), UncExpr([rhs],[ 1.],lhs.constant))
(-)(lhs::AffExpr, rhs::Uncertain) = UncVarExpr(lhs.vars, map(UncExpr,lhs.coeffs), UncExpr([rhs],[-1.],lhs.constant))
(*)(lhs::AffExpr, rhs::Uncertain) = UncVarExpr(lhs.vars,         rhs.*lhs.coeffs , UncExpr([rhs],[lhs.constant],0.0))
(/)(lhs::AffExpr, rhs::Uncertain) = error("Cannot divide $AFF by $UNC")
# AffExpr--UncExpr
(+)(lhs::AffExpr, rhs::UncExpr)  = UncVarExpr(lhs.vars, map(UncExpr,lhs.coeffs), lhs.constant+rhs)
(-)(lhs::AffExpr, rhs::UncExpr)  = UncVarExpr(lhs.vars, map(UncExpr,lhs.coeffs), lhs.constant-rhs)
(*)(lhs::AffExpr, rhs::UncExpr)  = UncVarExpr(lhs.vars,         rhs.*lhs.coeffs , lhs.constant*rhs)
(/)(lhs::AffExpr, rhs::UncExpr)  = error("Cannot divide $AFF by $UNE")
# AffExpr--UncVarExpr
(+)(lhs::AffExpr, rhs::UncVarExpr) = UncVarExpr(
  vcat(lhs.vars, rhs.vars),
  vcat(map(UncExpr,lhs.coeffs), rhs.coeffs),
  lhs.constant + rhs.constant)
(-)(lhs::AffExpr, rhs::UncVarExpr) = UncVarExpr(
  vcat(lhs.vars, rhs.vars),
  vcat(map(UncExpr,lhs.coeffs), -1.*rhs.coeffs),
  lhs.constant - rhs.constant)
(*)(lhs::AffExpr, rhs::UncVarExpr) = error("Cannot multiply $AFF by $UVE")
(/)(lhs::AffExpr, rhs::UncVarExpr) = error("Cannot divide $AFF by $UVE")
# AffExpr--GenericNormExpr{Uncertain}
(+){P}(lhs::AffExpr, rhs::GenericNormExpr{P,Float64,Uncertain}) =
    length(lhs.vars) == 0 ? lhs.constant + rhs : error("Cannot add $AFF by $UNM")
(-){P}(lhs::AffExpr, rhs::GenericNormExpr{P,Float64,Uncertain}) =
    length(lhs.vars) == 0 ? lhs.constant - rhs : error("Cannot substract $AFF by $UNM")
# AffExpr--Adaptive
(+)(a::AffExpr, x::Adaptive) = AdaptExpr(vcat(a.vars, x), vcat(a.coeffs,  1), a.constant)
(-)(a::AffExpr, x::Adaptive) = AdaptExpr(vcat(a.vars, x), vcat(a.coeffs, -1), a.constant)
(*)(a::AffExpr, x::Adaptive) = error("Cannot multiply $AFF by $ADP")
(/)(a::AffExpr, x::Adaptive) = error("Cannot divide $AFF by $ADP")
# AffExpr--AdaptiveAffExpr
(+)(a::AffExpr, b::AdaptExpr) = AdaptExpr(vcat(a.vars,   b.vars),
                                                vcat(a.coeffs, b.coeffs),
                                                a.constant + b.constant)
(-)(a::AffExpr, b::AdaptExpr) = AdaptExpr(vcat(a.vars,   b.vars),
                                                vcat(a.coeffs,-b.coeffs),
                                                a.constant - b.constant)
(*)(a::AffExpr, b::AdaptExpr) = error("Cannot multiply $AFF by $AFF")
(/)(a::AffExpr, b::AdaptExpr) = error("Cannot divide $AFF by $AFF")


#-----------------------------------------------------------------------
# 5. QuadExpr
# Nothing additional supported


#-----------------------------------------------------------------------
# 6. GenericNormExpr
# GenericNormExpr--Uncertain
(+){P,C,V<:Uncertain}(lhs::GenericNormExpr{P,C,V},rhs::Uncertain) = error("Cannot add $UNM to $UNC")
(-){P,C,V<:Uncertain}(lhs::GenericNormExpr{P,C,V},rhs::Uncertain) = error("Cannot subtract $UNM by $UNC")
(*){P,C,V<:Uncertain}(lhs::GenericNormExpr{P,C,V},rhs::Uncertain) = error("Cannot multiply $UNM by $UNC")
(/){P,C,V<:Uncertain}(lhs::GenericNormExpr{P,C,V},rhs::Uncertain) = error("Cannot divide $UNM by $UNC")
# GenericNormExpr--UncExpr
(+){P,C,V<:Uncertain}(lhs::GenericNormExpr{P,C,V},rhs::UncExpr) = error("Cannot add $UNM to $UNE")
(-){P,C,V<:Uncertain}(lhs::GenericNormExpr{P,C,V},rhs::UncExpr) = error("Cannot subtract $UNM by $UNE")
(*){P,C,V<:Uncertain}(lhs::GenericNormExpr{P,C,V},rhs::UncExpr) = error("Cannot multiply $UNM by $UNE")
(/){P,C,V<:Uncertain}(lhs::GenericNormExpr{P,C,V},rhs::UncExpr) = error("Cannot divide $UNM by $UNE")
# GenericNormExpr--UncVarExpr
(+){P,C,V<:Uncertain}(lhs::GenericNormExpr{P,C,V},rhs::UncVarExpr) = error("Cannot add $UNM to $UVE")
(-){P,C,V<:Uncertain}(lhs::GenericNormExpr{P,C,V},rhs::UncVarExpr) = error("Cannot subtract $UNM by $UVE")
(*){P,C,V<:Uncertain}(lhs::GenericNormExpr{P,C,V},rhs::UncVarExpr) = error("Cannot multiply $UNM by $UVE")
(/){P,C,V<:Uncertain}(lhs::GenericNormExpr{P,C,V},rhs::UncVarExpr) = error("Cannot divide $UNM by $UVE")


#-----------------------------------------------------------------------
# 7. Uncertain
(-)(lhs::Uncertain) = UncExpr(-1,lhs)
# Uncertain--Number
(+)(lhs::Uncertain, rhs::Number) = UncExpr([lhs],[1.0], rhs)
(-)(lhs::Uncertain, rhs::Number) = UncExpr([lhs],[1.0],-rhs)
(*)(lhs::Uncertain, rhs::Number) = UncExpr(  rhs, lhs)
(/)(lhs::Uncertain, rhs::Number) = UncExpr(1/rhs, lhs)
# Uncertain--Variable
(+)(lhs::Uncertain, rhs::Variable) = (+)(rhs, lhs)
(-)(lhs::Uncertain, rhs::Variable) = UncVarExpr([rhs],[UncExpr(-1)],UncExpr(lhs))
(*)(lhs::Uncertain, rhs::Variable) = (*)(rhs, lhs)
(/)(lhs::Uncertain, rhs::Variable) = error("Cannot divide $UNC by a variable")
# Uncertain--GenericNorm
(*){P,C,V<:Uncertain}(lhs::Uncertain, rhs::GenericNorm{P,C,V}) = error("Cannot multiply $UNC by $UNM")
(/){P,C,V<:Uncertain}(lhs::Uncertain, rhs::GenericNorm{P,C,V}) = error("Cannot multiply $UNC by $UNM")
# Uncertain--AffExpr
(+)(lhs::Uncertain, rhs::AffExpr) = (+)(rhs, lhs)
(-)(lhs::Uncertain, rhs::AffExpr) = UncVarExpr(rhs.vars, map(UncExpr,-rhs.coeffs), UncExpr([lhs],[1.0],-rhs.constant))
(*)(lhs::Uncertain, rhs::AffExpr) = (*)(rhs, lhs)
(/)(lhs::Uncertain, rhs::AffExpr) = error("Cannot divide $UNC by $AFF")
# Uncertain--GenericNormExpr
(+){P,C,V<:Uncertain}(lhs::Uncertain, rhs::GenericNormExpr{P,C,V}) = error("Cannot add $UNC to $UNM")
(-){P,C,V<:Uncertain}(lhs::Uncertain, rhs::GenericNormExpr{P,C,V}) = error("Cannot subtract $UNC by $UNM")
(*){P,C,V<:Uncertain}(lhs::Uncertain, rhs::GenericNormExpr{P,C,V}) = error("Cannot multiply $UNC by $UNM")
(/){P,C,V<:Uncertain}(lhs::Uncertain, rhs::GenericNormExpr{P,C,V}) = error("Cannot divide $UNC by $UNM")
# Uncertain--Uncertain
(+)(lhs::Uncertain, rhs::Uncertain) = UncExpr([lhs,rhs], Float64[1, 1], 0.0)
(-)(lhs::Uncertain, rhs::Uncertain) = UncExpr([lhs,rhs], Float64[1,-1], 0.0)
(*)(lhs::Uncertain, rhs::Uncertain) = error("Cannot multiply $UNC by $UNC")
(/)(lhs::Uncertain, rhs::Uncertain) = error("Cannot divide $UNC by $UNC")
# Uncertain--UncExpr
(+)(lhs::Uncertain, rhs::UncExpr) = UncExpr(vcat(lhs,rhs.vars),vcat(1.0, rhs.coeffs), rhs.constant)
(-)(lhs::Uncertain, rhs::UncExpr) = UncExpr(vcat(lhs,rhs.vars),vcat(1.0,-rhs.coeffs),-rhs.constant)
(*)(lhs::Uncertain, rhs::UncExpr) = error("Cannot multiply $UNC by $UNE")
(/)(lhs::Uncertain, rhs::UncExpr) = error("Cannot divide $UNC by $UNE")
# Uncertain--UncVarExpr
(+)(lhs::Uncertain, rhs::UncVarExpr) = UncVarExpr(rhs.vars,    rhs.coeffs,lhs+rhs.constant)
(-)(lhs::Uncertain, rhs::UncVarExpr) = UncVarExpr(rhs.vars,-1.*rhs.coeffs,lhs-rhs.constant)
(*)(lhs::Uncertain, rhs::UncVarExpr) = error("Cannot multiply $UNC by $UVE")
(/)(lhs::Uncertain, rhs::UncVarExpr) = error("Cannot divide $UNC by $UVE")
# Uncertain--Adaptive
(+)(u::Uncertain, x::Adaptive) = UncVarExpr([x], [ 1], u)
(-)(u::Uncertain, x::Adaptive) = UncVarExpr([x], [-1], u)
(*)(u::Uncertain, x::Adaptive) = UncVarExpr([x], [ u], 0)
(/)(u::Uncertain, x::Adaptive) = error("Cannot divide $UNC by $ADP")
# Uncertain--AdaptExpr
(+)(u::Uncertain, x::AdaptExpr) = UncVarExpr(copy(x.vars), copy(x.coeffs), u + x.constant)
(-)(u::Uncertain, x::AdaptExpr) = UncVarExpr(copy(x.vars),     -x.coeffs , u - x.constant)
(*)(u::Uncertain, x::AdaptExpr) = UncVarExpr(copy(x.vars), u .* x.coeffs , u * x.constant)
(/)(u::Uncertain, x::AdaptExpr) = error("Cannout divide $UNC by $AFF")


#-----------------------------------------------------------------------
# 8. UncExpr
# UncExpr--Number        - handled by JuMP
# UncExpr--Variable
(+)(lhs::UncExpr, rhs::Variable) = (+)(rhs,lhs)
(-)(lhs::UncExpr, rhs::Variable) = UncVarExpr([rhs],[UncExpr(-1)],lhs)
(*)(lhs::UncExpr, rhs::Variable) = (*)(rhs,lhs)
(/)(lhs::UncExpr, rhs::Variable) = error("Cannot divide $UNE by a variable")
# UncExpr--AffExpr
(+)(lhs::UncExpr, rhs::AffExpr) = (+)( rhs,lhs)
(-)(lhs::UncExpr, rhs::AffExpr) = (+)(-rhs,lhs)
(*)(lhs::UncExpr, rhs::AffExpr) = (*)( rhs,lhs)
(/)(lhs::UncExpr, rhs::AffExpr) = error("Cannot divide $UNE by $AFF")
# UncExpr--GenericNormExpr
(+){P,C,V<:Uncertain}(lhs::UncExpr, rhs::GenericNormExpr{P,C,V}) = error("Cannot add $UNE to $UNM")
(-){P,C,V<:Uncertain}(lhs::UncExpr, rhs::GenericNormExpr{P,C,V}) = error("Cannot subtract $UNE by $UNM")
(*){P,C,V<:Uncertain}(lhs::UncExpr, rhs::GenericNormExpr{P,C,V}) = error("Cannot multiply $UNE by $UNM")
(/){P,C,V<:Uncertain}(lhs::UncExpr, rhs::GenericNormExpr{P,C,V}) = error("Cannot divide $UNE by $UNM")
# UncExpr--Uncertain
(+)(lhs::UncExpr, rhs::Uncertain) = (+)(rhs,lhs)
(-)(lhs::UncExpr, rhs::Uncertain) = UncExpr(vcat(rhs,lhs.vars),vcat(-1.0,lhs.coeffs),lhs.constant)
(*)(lhs::UncExpr, rhs::Uncertain) = (*)(rhs,lhs)
(/)(lhs::UncExpr, rhs::Uncertain) = error("Cannot divide $UNE by $UNC")
# UncExpr--UncExpr
(*)(lhs::UncExpr, rhs::UncExpr) = error("Cannot multiply $UNE by $UNE")
(/)(lhs::UncExpr, rhs::UncExpr) = error("Cannot divide $UNE by $UNE")
# UncExpr--UncVarExpr
(+)(lhs::UncExpr, rhs::UncVarExpr) = UncVarExpr(rhs.vars,    rhs.coeffs,lhs+rhs.constant)
(-)(lhs::UncExpr, rhs::UncVarExpr) = UncVarExpr(rhs.vars,-1.*rhs.coeffs,lhs-rhs.constant)
(*)(lhs::UncExpr, rhs::UncVarExpr) = (length(lhs.vars) == 0) ?
                                        lhs.constant * rhs : # LHS is just a constant, so OK
                                        error("Cannot multiply $UNE by $UVE")
(/)(lhs::UncExpr, rhs::UncVarExpr) = error("Cannot divide $UNE by $UVE")
# UncExpr--Adaptive
(+)(u::UncExpr, x::Adaptive) = UncVarExpr([x], [ 1], u)
(-)(u::UncExpr, x::Adaptive) = UncVarExpr([x], [-1], u)
(*)(u::UncExpr, x::Adaptive) = UncVarExpr([x], [ u], 0)
(/)(u::UncExpr, x::Adaptive) = error("Cannot divide $UNE by $ADP")
# UncExpr--AdaptExpr
(+)(u::UncExpr, x::AdaptExpr) = UncVarExpr(copy(x.vars), copy(x.coeffs), u + x.constant)
(-)(u::UncExpr, x::AdaptExpr) = UncVarExpr(copy(x.vars),     -x.coeffs , u - x.constant)
(*)(u::UncExpr, x::AdaptExpr) = UncVarExpr(copy(x.vars), u .* x.coeffs , u * x.constant)
(/)(u::UncExpr, x::AdaptExpr) = error("Cannot divide $UNE by $AFF")


#-----------------------------------------------------------------------
# 9. UncVarExpr
# UncVarExpr--Number     - handled by JuMP
# UncVarExpr--Variable
(+)(lhs::UncVarExpr, rhs::Variable) = (+)(rhs,lhs)
(-)(lhs::UncVarExpr, rhs::Variable) = UncVarExpr(vcat(lhs.vars,rhs),vcat(lhs.coeffs,UncExpr(-1)), lhs.constant)
(*)(lhs::UncVarExpr, rhs::Variable) = error("Cannot multiply $UVE by a variable")
(/)(lhs::UncVarExpr, rhs::Variable) = error("Cannot divide $UVE by a variable")
# UncVarExpr--AffExpr
(+)(lhs::UncVarExpr, rhs::AffExpr) = (+)(rhs,lhs)
(-)(lhs::UncVarExpr, rhs::AffExpr) = UncVarExpr(
  vcat(lhs.vars,    rhs.vars),
  vcat(lhs.coeffs, -rhs.coeffs),
  lhs.constant - rhs.constant)
(*)(lhs::UncVarExpr, rhs::AffExpr) = error("Cannot multiply $UVE by $AFF")
(/)(lhs::UncVarExpr, rhs::AffExpr) = error("Cannot divide $UVE by $AFF")
# UncVarExpr--GenericNormExpr
(+){P,C,V<:Uncertain}(lhs::UncVarExpr, rhs::GenericNormExpr{P,C,V}) = error("Cannot add $UVE to $UNM")
(-){P,C,V<:Uncertain}(lhs::UncVarExpr, rhs::GenericNormExpr{P,C,V}) = error("Cannot subtract $UVE by $UNM")
(*){P,C,V<:Uncertain}(lhs::UncVarExpr, rhs::GenericNormExpr{P,C,V}) = error("Cannot multiply $UVE by $UNM")
(/){P,C,V<:Uncertain}(lhs::UncVarExpr, rhs::GenericNormExpr{P,C,V}) = error("Cannot divide $UVE by $UNM")
# UncVarExpr--Uncertain
(+)(lhs::UncVarExpr, rhs::Uncertain) = UncVarExpr(lhs.vars,lhs.coeffs,lhs.constant+rhs)
(-)(lhs::UncVarExpr, rhs::Uncertain) = UncVarExpr(lhs.vars,lhs.coeffs,lhs.constant-rhs)
(*)(lhs::UncVarExpr, rhs::Uncertain) = error("Cannot multiply $UVE by $UNC")
(/)(lhs::UncVarExpr, rhs::Uncertain) = error("Cannot divide $UVE by $UNC")
# UncVarExpr--UncExpr
(+)(lhs::UncVarExpr, rhs::UncExpr) = UncVarExpr(lhs.vars,lhs.coeffs,lhs.constant+rhs)
(-)(lhs::UncVarExpr, rhs::UncExpr) = UncVarExpr(lhs.vars,lhs.coeffs,lhs.constant-rhs)
(*)(lhs::UncVarExpr, rhs::UncExpr) = error("Cannot multiply $UVE by $UNE")
(/)(lhs::UncVarExpr, rhs::UncExpr) = error("Cannot divide $UVE by $UNE")
# UncVarExpr--UncVarExpr
(*)(lhs::UncVarExpr, rhs::UncVarExpr) = error("Cannot multiply $UVE by $UVE")
(/)(lhs::UncVarExpr, rhs::UncVarExpr) = error("Cannot divide $UVE by $UVE")
# UncVarExpr--Adaptive
+(a::UncVarExpr, x::Adaptive) = UncVarExpr(vcat(a.vars, x),
                                           vcat(a.coeffs, UncExpr(1)),
                                           a.constant)
-(a::UncVarExpr, x::Adaptive) = UncVarExpr(vcat(a.vars, x),
                                           vcat(a.coeffs, UncExpr(-1)),
                                           a.constant)
*(a::UncVarExpr, x::Adaptive) = error("Cannot multiply $UVE by $ADP")
/(a::UncVarExpr, x::Adaptive) = error("Cannot divide $UVE by $ADP")
# UncVarExpr--AdaptExpr
+(a::UncVarExpr, b::AdaptExpr) = UncVarExpr(vcat(a.vars, b.vars),
                                               vcat(a.coeffs, map(UncExpr, b.coeffs)),
                                               a.constant + b.constant)
-(a::UncVarExpr, b::AdaptExpr) = UncVarExpr(vcat(a.vars, b.vars),
                                               vcat(a.coeffs, map(UncExpr,-b.coeffs)),
                                               a.constant - b.constant)
*(a::UncVarExpr, b::AdaptExpr) = error("Cannot multiply $UVE by $AFF")
/(a::UncVarExpr, b::AdaptExpr) = error("Cannot divide $UVE by $AFF")


#-----------------------------------------------------------------------
# 10. Adaptive
# Adaptive--Number
(+)(x::Adaptive, c::Number) = +(  c, x)
(-)(x::Adaptive, c::Number) = +( -c, x)
(*)(x::Adaptive, c::Number) = *(  c, x)
(/)(x::Adaptive, c::Number) = *(1/c, x)
# Adaptive--Variable
(+)(x::Adaptive, v::Variable) = AdaptExpr([x,v], [1,+1], 0)
(-)(x::Adaptive, v::Variable) = AdaptExpr([x,v], [1,-1], 0)
(*)(x::Adaptive, v::Variable) = error("Cannot multiply $ADP by a variable")
(/)(x::Adaptive, v::Variable) = error("Cannot divide $ADP by a variable")
# Adaptive--AffExpr
(+)(x::Adaptive, a::AffExpr) = AdaptExpr(vcat(x, a.vars), vcat(1,  a.coeffs),  a.constant)
(-)(x::Adaptive, a::AffExpr) = AdaptExpr(vcat(x, a.vars), vcat(1, -a.coeffs), -a.constant)
(*)(x::Adaptive, a::AffExpr) = error("Cannot multiply $ADP by $AFF")
(/)(x::Adaptive, a::AffExpr) = error("Cannot divide $ADP by $AFF")
# Adaptive--Uncertain
(+)(x::Adaptive, u::Uncertain) = UncVarExpr([x], [1],  u)
(-)(x::Adaptive, u::Uncertain) = UncVarExpr([x], [1], -u)
(*)(x::Adaptive, u::Uncertain) = UncVarExpr([x], [u],  0)
(/)(x::Adaptive, u::Uncertain) = error("Cannot divide $ADP by $UNC")
# Adaptive--UncExpr
(+)(x::Adaptive, u::UncExpr) = UncVarExpr([x], [1],  u)
(-)(x::Adaptive, u::UncExpr) = UncVarExpr([x], [1], -u)
(*)(x::Adaptive, u::UncExpr) = UncVarExpr([x], [u],  0)
(/)(x::Adaptive, u::UncExpr) = error("Cannot divide $ADP by $UNE")
# Adaptive--UncVarExpr
(+)(x::Adaptive, u::UncVarExpr) = UncVarExpr(vcat(x, u.vars), vcat(1,  u.coeffs),  u.constant)
(-)(x::Adaptive, u::UncVarExpr) = UncVarExpr(vcat(x, u.vars), vcat(1, -u.coeffs), -u.constant)
(*)(x::Adaptive, u::UncVarExpr) = error("Cannot multiply $ADP by $UVE")
(/)(x::Adaptive, u::UncVarExpr) = error("Cannot divide $ADP by $UVE")
# Adaptive--Adaptive
(+)(a::Adaptive, b::Adaptive) = AdaptExpr(Adaptive[a,b], Float64[1, 1], 0)
(-)(a::Adaptive, b::Adaptive) = AdaptExpr(Adaptive[a,b], Float64[1,-1], 0)
(*)(a::Adaptive, b::Adaptive) = error("Cannot multiply $ADP by $ADP")
(/)(a::Adaptive, b::Adaptive) = error("Cannot divide $ADP by $ADP")
# Adaptive--AdaptExpr
(+)(a::Adaptive, b::AdaptExpr) = AdaptExpr(vcat(a, b.vars), vcat(1,  b.coeffs),  b.constant)
(-)(a::Adaptive, b::AdaptExpr) = AdaptExpr(vcat(a, b.vars), vcat(1, -b.coeffs), -b.constant)
(*)(a::Adaptive, b::AdaptExpr) = error("Cannot multiply $ADP by $AFF")
(/)(a::Adaptive, b::AdaptExpr) = error("Cannot divide $ADP by $AFF")


#-----------------------------------------------------------------------
# 11. AdaptExpr
# AdaptExpr--Number
(+)(a::AdaptExpr, c::Number) = AdaptExpr(copy(a.vars), copy(a.coeffs), a.constant + c)
(-)(a::AdaptExpr, c::Number) = AdaptExpr(copy(a.vars), copy(a.coeffs), a.constant - c)
(*)(a::AdaptExpr, c::Number) = AdaptExpr(copy(a.vars), copy(a.coeffs)*c, a.constant*c)
(/)(a::AdaptExpr, c::Number) = AdaptExpr(copy(a.vars), copy(a.coeffs)/c, a.constant/c)
# AdaptExpr--Variable
(+)(x::AdaptExpr, v::Variable) = AdaptExpr(vcat(x.vars, v), vcat(x.coeffs,  1), x.constant)
(-)(x::AdaptExpr, v::Variable) = AdaptExpr(vcat(x.vars, v), vcat(x.coeffs, -1), x.constant)
(*)(x::AdaptExpr, v::Variable) = error("Cannot multiply $AFF by a variable")
(/)(x::AdaptExpr, v::Variable) = error("Cannot divide $AFF by a variable")
# AdaptExpr--AffExpr
(+)(x::AdaptExpr, a::AffExpr) = AdaptExpr(vcat(x.vars, a.vars), vcat(x.coeffs,  a.coeffs), x.constant + a.constant)
(-)(x::AdaptExpr, a::AffExpr) = AdaptExpr(vcat(x.vars, a.vars), vcat(x.coeffs, -a.coeffs), x.constant - a.constant)
(*)(x::AdaptExpr, a::AffExpr) = error("Cannot multiply $AFF by $AFF")
(/)(x::AdaptExpr, a::AffExpr) = error("Cannot divide $AFF by $AFF")
# AdaptExpr--Uncertain
(+)(a::AdaptExpr, u::Uncertain) = UncVarExpr(copy(a.vars), copy(a.coeffs), a.constant + u)
(-)(a::AdaptExpr, u::Uncertain) = UncVarExpr(copy(a.vars), copy(a.coeffs), a.constant - u)
(*)(a::AdaptExpr, u::Uncertain) = UncVarExpr(copy(a.vars), copy(a.coeffs) * u, a.constant * u)
(/)(a::AdaptExpr, u::Uncertain) = error("Cannot divide $AFF by $UNC")
# AdaptExpr--UncExpr
(+)(a::AdaptExpr, b::UncExpr) = +( b, a)
(-)(a::AdaptExpr, b::UncExpr) = +(-b, a)
(*)(a::AdaptExpr, b::UncExpr) = *( b, a)
(/)(a::AdaptExpr, b::UncExpr) = /( b, a)  # not quite
# AdaptExpr--UncVarExpr
(+)(a::AdaptExpr, b::UncVarExpr) = UncVarExpr(vcat(a.vars, b.vars),
                                                 vcat(map(UncExpr,a.coeffs), b.coeffs),
                                                 a.constant + b.constant)
(-)(a::AdaptExpr, b::UncVarExpr) = UncVarExpr(vcat(a.vars, b.vars),
                                                 vcat(map(UncExpr,a.coeffs), -b.coeffs),
                                                 a.constant - b.constant)
(*)(a::AdaptExpr, b::UncVarExpr) = error("Cannot multiply $AFF by $UVE")
(/)(a::AdaptExpr, b::UncVarExpr) = error("Cannot divide $AFF by $UVE")
# AdaptExpr--Adaptive
(+)(a::AdaptExpr, b::Adaptive) = AdaptExpr(vcat(a.vars,b), vcat(a.coeffs, 1), a.constant)
(-)(a::AdaptExpr, b::Adaptive) = AdaptExpr(vcat(a.vars,b), vcat(a.coeffs,-1), a.constant)
(*)(a::AdaptExpr, b::Adaptive) = error("Cannot multiply $AFF by $ADP")
(/)(a::AdaptExpr, b::Adaptive) = error("Cannot divide $AFF by $ADP")
# AdaptExpr--AdaptExpr
(+)(a::AdaptExpr, b::AdaptExpr) = AdaptExpr(vcat(a, b.vars), vcat(1,  b.coeffs), a.constant+b.constant)
(-)(a::AdaptExpr, b::AdaptExpr) = AdaptExpr(vcat(a, b.vars), vcat(1, -b.coeffs), a.constant-b.constant)
(*)(a::AdaptExpr, b::AdaptExpr) = error("Cannot multiply $AFF by $AFF")
(/)(a::AdaptExpr, b::AdaptExpr) = error("Cannot divide $AFF by $AFF")


#-----------------------------------------------------------------------
# For matrix operations
Base.promote_rule{R<:Real}(::Type{R},         ::Type{Uncertain} ) = UncExpr
Base.promote_rule{R<:Real}(::Type{R},         ::Type{UncExpr}   ) = UncExpr
Base.promote_rule{R<:Real}(::Type{R},         ::Type{UncVarExpr}) = UncVarExpr
Base.promote_rule{R<:Real}(::Type{R},         ::Type{Adaptive})   = AdaptExpr
Base.promote_rule{R<:Real}(::Type{R},         ::Type{AdaptExpr})  = AdaptExpr

Base.promote_rule{R<:Real}(::Type{Uncertain}, ::Type{R})          = UncExpr
Base.promote_rule(         ::Type{Uncertain}, ::Type{Uncertain})  = UncExpr
Base.promote_rule(         ::Type{Uncertain}, ::Type{UncExpr})    = UncExpr
Base.promote_rule(         ::Type{Uncertain}, ::Type{UncVarExpr}) = UncVarExpr
Base.promote_rule(         ::Type{Uncertain}, ::Type{Adaptive})   = UncVarExpr
Base.promote_rule(         ::Type{Uncertain}, ::Type{AdaptExpr})  = UncVarExpr

Base.promote_rule{R<:Real}(::Type{UncExpr},   ::Type{R})          = UncExpr
Base.promote_rule(         ::Type{UncExpr},   ::Type{Uncertain})  = UncExpr
Base.promote_rule(         ::Type{UncExpr},   ::Type{UncExpr})    = UncExpr
Base.promote_rule(         ::Type{UncExpr},   ::Type{UncVarExpr}) = UncVarExpr
Base.promote_rule(         ::Type{UncExpr},   ::Type{Adaptive})   = UncVarExpr
Base.promote_rule(         ::Type{UncExpr},   ::Type{AdaptExpr})  = UncVarExpr

Base.promote_rule{R<:Real}(::Type{UncVarExpr},::Type{R})          = UncVarExpr
Base.promote_rule(         ::Type{UncVarExpr},::Type{Uncertain})  = UncVarExpr
Base.promote_rule(         ::Type{UncVarExpr},::Type{UncExpr})    = UncVarExpr
Base.promote_rule(         ::Type{UncVarExpr},::Type{UncVarExpr}) = UncVarExpr
Base.promote_rule(         ::Type{UncVarExpr},::Type{Adaptive})   = UncVarExpr
Base.promote_rule(         ::Type{UncVarExpr},::Type{AdaptExpr})  = UncVarExpr

Base.promote_rule{R<:Real}(::Type{Adaptive},  ::Type{R})          = AdaptExpr
Base.promote_rule(         ::Type{Adaptive},  ::Type{Uncertain})  = UncVarExpr
Base.promote_rule(         ::Type{Adaptive},  ::Type{UncExpr})    = UncVarExpr
Base.promote_rule(         ::Type{Adaptive},  ::Type{UncVarExpr}) = UncVarExpr
Base.promote_rule(         ::Type{Adaptive},  ::Type{Adaptive})   = AdaptExpr
Base.promote_rule(         ::Type{Adaptive},  ::Type{AdaptExpr})  = AdaptExpr

Base.promote_rule{R<:Real}(::Type{AdaptExpr}, ::Type{R})          = AdaptExpr
Base.promote_rule(         ::Type{AdaptExpr}, ::Type{Uncertain})  = UncVarExpr
Base.promote_rule(         ::Type{AdaptExpr}, ::Type{UncExpr})    = UncVarExpr
Base.promote_rule(         ::Type{AdaptExpr}, ::Type{UncVarExpr}) = UncVarExpr
Base.promote_rule(         ::Type{AdaptExpr}, ::Type{Adaptive})   = AdaptExpr
Base.promote_rule(         ::Type{AdaptExpr}, ::Type{AdaptExpr})  = AdaptExpr
