#-----------------------------------------------------------------------
# JuMPeR  --  JuMP Extension for Robust Optimization
# http://github.com/IainNZ/JuMPeR.jl
#-----------------------------------------------------------------------
# Copyright (c) 2015: Iain Dunning
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
# 4. [Generic]AffExpr
# 5. QuadExpr <- We don't support any interactions with QuadExpr
# 6. [Generic]NormExpr
# 7. Uncertain
# 8. UncExpr
# 9. FullAffExpr
#-----------------------------------------------------------------------

# Define error message names for the sake of consistency
const AFF = "an affine function of variables"
const UNC = "an uncertain parameter"
const UAE = "an affine function of uncertain parameters"
const FAE = "an affine function of variables and uncertain parameters"
const UNM = "a norm of uncertain parameters"

#-----------------------------------------------------------------------
# 1. Number
# Number--Uncertain
(+)(lhs::Number, rhs::Uncertain) = UncExpr([rhs], [  1], lhs)
(-)(lhs::Number, rhs::Uncertain) = UncExpr([rhs], [ -1], lhs)
(*)(lhs::Number, rhs::Uncertain) = UncExpr(lhs, rhs)
(/)(lhs::Number, rhs::Uncertain) = error("Cannot divide a number by $UNC")
# Number--UncExpr        - handled by JuMP
# Number--FullAffExpr     - handled by JuMP

#-----------------------------------------------------------------------
# 2. Variable
# Variable--Uncertain
(+)(lhs::Variable, rhs::Uncertain) = FullAffExpr([lhs],[UncExpr(  1)], UncExpr(rhs))
(-)(lhs::Variable, rhs::Uncertain) = FullAffExpr([lhs],[UncExpr(  1)],-UncExpr(rhs))
(*)(lhs::Variable, rhs::Uncertain) = FullAffExpr([lhs],[UncExpr(rhs)], UncExpr())
(/)(lhs::Variable, rhs::Uncertain) = error("Cannot divide a variable by $UNC")
# Variable--UncExpr
(+)(lhs::Variable, rhs::UncExpr) = FullAffExpr([lhs],[UncExpr(1)],       rhs)
(-)(lhs::Variable, rhs::UncExpr) = FullAffExpr([lhs],[UncExpr(1)],      -rhs)
(*)(lhs::Variable, rhs::UncExpr) = FullAffExpr([lhs],[        rhs],UncExpr())
(/)(lhs::Variable, rhs::UncExpr) = error("Cannot divide a variable by $UAE")
# Variable--FullAffExpr
(+)(lhs::Variable, rhs::FullAffExpr) = FullAffExpr(vcat(rhs.vars,lhs),vcat(    rhs.coeffs,UncExpr(1)), rhs.constant)
(-)(lhs::Variable, rhs::FullAffExpr) = FullAffExpr(vcat(rhs.vars,lhs),vcat(-1.*rhs.coeffs,UncExpr(1)),-rhs.constant)
(*)(lhs::Variable, rhs::FullAffExpr) = error("Cannot multiply a variable by $FAE")
(/)(lhs::Variable, rhs::FullAffExpr) = error("Cannot divide a variable by $FAE")

#-----------------------------------------------------------------------
# 3. Norm


#-----------------------------------------------------------------------
# 4. AffExpr
# AffExpr--Uncertain
(+)(lhs::AffExpr, rhs::Uncertain) = FullAffExpr(lhs.vars, map(UncExpr,lhs.coeffs), UncExpr([rhs],[ 1.],lhs.constant))
(-)(lhs::AffExpr, rhs::Uncertain) = FullAffExpr(lhs.vars, map(UncExpr,lhs.coeffs), UncExpr([rhs],[-1.],lhs.constant))
(*)(lhs::AffExpr, rhs::Uncertain) = FullAffExpr(lhs.vars,         rhs.*lhs.coeffs , UncExpr([rhs],[lhs.constant],0.0))
(/)(lhs::AffExpr, rhs::Uncertain) = error("Cannot divide $AFF by $UNC")
# AffExpr-UncExpr
(+)(lhs::AffExpr, rhs::UncExpr)  = FullAffExpr(lhs.vars, map(UncExpr,lhs.coeffs), lhs.constant+rhs)
(-)(lhs::AffExpr, rhs::UncExpr)  = FullAffExpr(lhs.vars, map(UncExpr,lhs.coeffs), lhs.constant-rhs)
(*)(lhs::AffExpr, rhs::UncExpr)  = FullAffExpr(lhs.vars,         rhs.*lhs.coeffs , lhs.constant*rhs)
(/)(lhs::AffExpr, rhs::UncExpr)  = error("Cannot divide $AFF by $UAE")
# AffExpr-FullAffExpr
(+)(lhs::AffExpr, rhs::FullAffExpr) = FullAffExpr(
  vcat(lhs.vars, rhs.vars),
  vcat(map(UncExpr,lhs.coeffs), rhs.coeffs),
  lhs.constant + rhs.constant)
(-)(lhs::AffExpr, rhs::FullAffExpr) = FullAffExpr(
  vcat(lhs.vars, rhs.vars),
  vcat(map(UncExpr,lhs.coeffs), -1.*rhs.coeffs),
  lhs.constant - rhs.constant)
(*)(lhs::AffExpr, rhs::FullAffExpr) = error("Cannot multiply $AFF by $FAE")
(/)(lhs::AffExpr, rhs::FullAffExpr) = error("Cannot divide $AFF by $FAE")
# AffExpr -- GenericNormExpr{Uncertain}
(+){P}(lhs::AffExpr, rhs::GenericNormExpr{P,Float64,Uncertain}) =
    length(lhs.vars) == 0 ? lhs.constant + rhs : error("Cannot add $AFF by $UNM")
(-){P}(lhs::AffExpr, rhs::GenericNormExpr{P,Float64,Uncertain}) =
    length(lhs.vars) == 0 ? lhs.constant - rhs : error("Cannot substract $AFF by $UNM")


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
(+){P,C,V<:Uncertain}(lhs::GenericNormExpr{P,C,V},rhs::UncExpr) = error("Cannot add $UNM to $UAE")
(-){P,C,V<:Uncertain}(lhs::GenericNormExpr{P,C,V},rhs::UncExpr) = error("Cannot subtract $UNM by $UAE")
(*){P,C,V<:Uncertain}(lhs::GenericNormExpr{P,C,V},rhs::UncExpr) = error("Cannot multiply $UNM by $UAE")
(/){P,C,V<:Uncertain}(lhs::GenericNormExpr{P,C,V},rhs::UncExpr) = error("Cannot divide $UNM by $UAE")
# GenericNormExpr--FullAffExpr
(+){P,C,V<:Uncertain}(lhs::GenericNormExpr{P,C,V},rhs::FullAffExpr) = error("Cannot add $UNM to $FAE")
(-){P,C,V<:Uncertain}(lhs::GenericNormExpr{P,C,V},rhs::FullAffExpr) = error("Cannot subtract $UNM by $FAE")
(*){P,C,V<:Uncertain}(lhs::GenericNormExpr{P,C,V},rhs::FullAffExpr) = error("Cannot multiply $UNM by $FAE")
(/){P,C,V<:Uncertain}(lhs::GenericNormExpr{P,C,V},rhs::FullAffExpr) = error("Cannot divide $UNM by $FAE")


#-----------------------------------------------------------------------
# Uncertain
(-)(lhs::Uncertain) = UncExpr(-1,lhs)
# Uncertain--Number
(+)(lhs::Uncertain, rhs::Number) = UncExpr([lhs],[1.0], rhs)
(-)(lhs::Uncertain, rhs::Number) = UncExpr([lhs],[1.0],-rhs)
(*)(lhs::Uncertain, rhs::Number) = UncExpr(  rhs, lhs)
(/)(lhs::Uncertain, rhs::Number) = UncExpr(1/rhs, lhs)
# Uncertain--Variable
(+)(lhs::Uncertain, rhs::Variable) = (+)(rhs, lhs)
(-)(lhs::Uncertain, rhs::Variable) = FullAffExpr([rhs],[UncExpr(-1)],UncExpr(lhs))
(*)(lhs::Uncertain, rhs::Variable) = (*)(rhs, lhs)
(/)(lhs::Uncertain, rhs::Variable) = error("Cannot divide $UNC by a variable")
# Uncertain--AffExpr
(+)(lhs::Uncertain, rhs::AffExpr) = (+)(rhs, lhs)
(-)(lhs::Uncertain, rhs::AffExpr) = FullAffExpr(rhs.vars, map(UncExpr,-rhs.coeffs), UncExpr([lhs],[1.0],-rhs.constant))
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
(*)(lhs::Uncertain, rhs::UncExpr) = error("Cannot multiply $UNC by $UAE")
(/)(lhs::Uncertain, rhs::UncExpr) = error("Cannot divide $UNC by $UAE")
# Uncertain--FullAffExpr
(+)(lhs::Uncertain, rhs::FullAffExpr) = FullAffExpr(rhs.vars,    rhs.coeffs,lhs+rhs.constant)
(-)(lhs::Uncertain, rhs::FullAffExpr) = FullAffExpr(rhs.vars,-1.*rhs.coeffs,lhs-rhs.constant)
(*)(lhs::Uncertain, rhs::FullAffExpr) = error("Cannot multiply $UNC by $FAE")
(/)(lhs::Uncertain, rhs::FullAffExpr) = error("Cannot divide $UNC by $FAE")

# UncExpr
# UncExpr--Number        - handled by JuMP
# UncExpr--Variable
(+)(lhs::UncExpr, rhs::Variable) = (+)(rhs,lhs)
(-)(lhs::UncExpr, rhs::Variable) = FullAffExpr([rhs],[UncExpr(-1)],lhs)
(*)(lhs::UncExpr, rhs::Variable) = (*)(rhs,lhs)
(/)(lhs::UncExpr, rhs::Variable) = error("Cannot divide $UAE by a variable")
# UncExpr--AffExpr
(+)(lhs::UncExpr, rhs::AffExpr) = (+)( rhs,lhs)
(-)(lhs::UncExpr, rhs::AffExpr) = (+)(-rhs,lhs)
(*)(lhs::UncExpr, rhs::AffExpr) = (*)( rhs,lhs)
(/)(lhs::UncExpr, rhs::AffExpr) = error("Cannot divide $UAE by $AFF")
# UncExpr--GenericNormExpr
(+){P,C,V<:Uncertain}(lhs::UncExpr, rhs::GenericNormExpr{P,C,V}) = error("Cannot add $UAE to $UNM")
(-){P,C,V<:Uncertain}(lhs::UncExpr, rhs::GenericNormExpr{P,C,V}) = error("Cannot subtract $UAE by $UNM")
(*){P,C,V<:Uncertain}(lhs::UncExpr, rhs::GenericNormExpr{P,C,V}) = error("Cannot multiply $UAE by $UNM")
(/){P,C,V<:Uncertain}(lhs::UncExpr, rhs::GenericNormExpr{P,C,V}) = error("Cannot divide $UAE by $UNM")
# UncExpr--Uncertain
(+)(lhs::UncExpr, rhs::Uncertain) = (+)(rhs,lhs)
(-)(lhs::UncExpr, rhs::Uncertain) = UncExpr(vcat(rhs,lhs.vars),vcat(-1.0,lhs.coeffs),lhs.constant)
(*)(lhs::UncExpr, rhs::Uncertain) = (*)(rhs,lhs)
(/)(lhs::UncExpr, rhs::Uncertain) = error("Cannot divide $UAE by $UNC")
# UncExpr--UncExpr
(*)(lhs::UncExpr, rhs::UncExpr) = error("Cannot multiply $UAE by $UAE")
(/)(lhs::UncExpr, rhs::UncExpr) = error("Cannot divide $UAE by $UAE")
# UncExpr--FullAffExpr
(+)(lhs::UncExpr, rhs::FullAffExpr) = FullAffExpr(rhs.vars,    rhs.coeffs,lhs+rhs.constant)
(-)(lhs::UncExpr, rhs::FullAffExpr) = FullAffExpr(rhs.vars,-1.*rhs.coeffs,lhs-rhs.constant)
(*)(lhs::UncExpr, rhs::FullAffExpr) = error("Cannot multiply $UAE by $FAE")
(/)(lhs::UncExpr, rhs::FullAffExpr) = error("Cannot divide $UAE by $FAE")

# FullAffExpr
# FullAffExpr--Number     - handled by JuMP
# FullAffExpr--Variable
(+)(lhs::FullAffExpr, rhs::Variable) = (+)(rhs,lhs)
(-)(lhs::FullAffExpr, rhs::Variable) = FullAffExpr(vcat(lhs.vars,rhs),vcat(lhs.coeffs,UncExpr(-1)), lhs.constant)
(*)(lhs::FullAffExpr, rhs::Variable) = error("Cannot multiply $FAE by a variable")
(/)(lhs::FullAffExpr, rhs::Variable) = error("Cannot divide $FAE by a variable")
# FullAffExpr--AffExpr
(+)(lhs::FullAffExpr, rhs::AffExpr) = (+)(rhs,lhs)
(-)(lhs::FullAffExpr, rhs::AffExpr) = FullAffExpr(
  vcat(lhs.vars,    rhs.vars),
  vcat(lhs.coeffs, -rhs.coeffs),
  lhs.constant - rhs.constant)
(*)(lhs::FullAffExpr, rhs::AffExpr) = error("Cannot multiply $FAE by $AFF")
(/)(lhs::FullAffExpr, rhs::AffExpr) = error("Cannot divide $FAE by $AFF")
# FullAffExpr--GenericNormExpr
(+){P,C,V<:Uncertain}(lhs::FullAffExpr, rhs::GenericNormExpr{P,C,V}) = error("Cannot add $FAE to $UNM")
(-){P,C,V<:Uncertain}(lhs::FullAffExpr, rhs::GenericNormExpr{P,C,V}) = error("Cannot subtract $FAE by $UNM")
(*){P,C,V<:Uncertain}(lhs::FullAffExpr, rhs::GenericNormExpr{P,C,V}) = error("Cannot multiply $FAE by $UNM")
(/){P,C,V<:Uncertain}(lhs::FullAffExpr, rhs::GenericNormExpr{P,C,V}) = error("Cannot divide $FAE by $UNM")
# FullAffExpr--Uncertain
(+)(lhs::FullAffExpr, rhs::Uncertain) = FullAffExpr(lhs.vars,lhs.coeffs,lhs.constant+rhs)
(-)(lhs::FullAffExpr, rhs::Uncertain) = FullAffExpr(lhs.vars,lhs.coeffs,lhs.constant-rhs)
(*)(lhs::FullAffExpr, rhs::Uncertain) = error("Cannot multiply $FAE by $UNC")
(/)(lhs::FullAffExpr, rhs::Uncertain) = error("Cannot divide $FAE by $UNC")
# FullAffExpr--UncExpr
(+)(lhs::FullAffExpr, rhs::UncExpr) = FullAffExpr(lhs.vars,lhs.coeffs,lhs.constant+rhs)
(-)(lhs::FullAffExpr, rhs::UncExpr) = FullAffExpr(lhs.vars,lhs.coeffs,lhs.constant-rhs)
(*)(lhs::FullAffExpr, rhs::UncExpr) = error("Cannot multiply $FAE by $UAE")
(/)(lhs::FullAffExpr, rhs::UncExpr) = error("Cannot divide $FAE by $UAE")
# FullAffExpr--FullAffExpr
(*)(lhs::FullAffExpr, rhs::FullAffExpr) = error("Cannot multiply $FAE by $FAE")
(/)(lhs::FullAffExpr, rhs::FullAffExpr) = error("Cannot divide $FAE by $FAE")


#-----------------------------------------------------------------------
# High-level operators like sum and dot are handled by JuMP
if VERSION < v"0.4-"
    import JuMP.vecdot
    export vecdot
end