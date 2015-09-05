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
# 9. UncAffExpr
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
# Number--UncAffExpr     - handled by JuMP

#-----------------------------------------------------------------------
# 2. Variable
# Variable--Uncertain
(+)(lhs::Variable, rhs::Uncertain) = UncAffExpr([lhs],[UncExpr(  1)], UncExpr(rhs))
(-)(lhs::Variable, rhs::Uncertain) = UncAffExpr([lhs],[UncExpr(  1)],-UncExpr(rhs))
(*)(lhs::Variable, rhs::Uncertain) = UncAffExpr([lhs],[UncExpr(rhs)], UncExpr())
(/)(lhs::Variable, rhs::Uncertain) = error("Cannot divide a variable by $UNC")
# Variable--UncExpr
(+)(lhs::Variable, rhs::UncExpr) = UncAffExpr([lhs],[UncExpr(1)],       rhs)
(-)(lhs::Variable, rhs::UncExpr) = UncAffExpr([lhs],[UncExpr(1)],      -rhs)
(*)(lhs::Variable, rhs::UncExpr) = UncAffExpr([lhs],[        rhs],UncExpr())
(/)(lhs::Variable, rhs::UncExpr) = error("Cannot divide a variable by $UAE")
# Variable--UncAffExpr
(+)(lhs::Variable, rhs::UncAffExpr) = UncAffExpr(vcat(rhs.vars,lhs),vcat(    rhs.coeffs,UncExpr(1)), rhs.constant)
(-)(lhs::Variable, rhs::UncAffExpr) = UncAffExpr(vcat(rhs.vars,lhs),vcat(-1.*rhs.coeffs,UncExpr(1)),-rhs.constant)
(*)(lhs::Variable, rhs::UncAffExpr) = error("Cannot multiply a variable by $FAE")
(/)(lhs::Variable, rhs::UncAffExpr) = error("Cannot divide a variable by $FAE")

#-----------------------------------------------------------------------
# 3. Norm


#-----------------------------------------------------------------------
# 4. AffExpr
# AffExpr--Uncertain
(+)(lhs::AffExpr, rhs::Uncertain) = UncAffExpr(lhs.vars, map(UncExpr,lhs.coeffs), UncExpr([rhs],[ 1.],lhs.constant))
(-)(lhs::AffExpr, rhs::Uncertain) = UncAffExpr(lhs.vars, map(UncExpr,lhs.coeffs), UncExpr([rhs],[-1.],lhs.constant))
(*)(lhs::AffExpr, rhs::Uncertain) = UncAffExpr(lhs.vars,         rhs.*lhs.coeffs , UncExpr([rhs],[lhs.constant],0.0))
(/)(lhs::AffExpr, rhs::Uncertain) = error("Cannot divide $AFF by $UNC")
# AffExpr-UncExpr
(+)(lhs::AffExpr, rhs::UncExpr)  = UncAffExpr(lhs.vars, map(UncExpr,lhs.coeffs), lhs.constant+rhs)
(-)(lhs::AffExpr, rhs::UncExpr)  = UncAffExpr(lhs.vars, map(UncExpr,lhs.coeffs), lhs.constant-rhs)
(*)(lhs::AffExpr, rhs::UncExpr)  = UncAffExpr(lhs.vars,         rhs.*lhs.coeffs , lhs.constant*rhs)
(/)(lhs::AffExpr, rhs::UncExpr)  = error("Cannot divide $AFF by $UAE")
# AffExpr-UncAffExpr
(+)(lhs::AffExpr, rhs::UncAffExpr) = UncAffExpr(
  vcat(lhs.vars, rhs.vars),
  vcat(map(UncExpr,lhs.coeffs), rhs.coeffs),
  lhs.constant + rhs.constant)
(-)(lhs::AffExpr, rhs::UncAffExpr) = UncAffExpr(
  vcat(lhs.vars, rhs.vars),
  vcat(map(UncExpr,lhs.coeffs), -1.*rhs.coeffs),
  lhs.constant - rhs.constant)
(*)(lhs::AffExpr, rhs::UncAffExpr) = error("Cannot multiply $AFF by $FAE")
(/)(lhs::AffExpr, rhs::UncAffExpr) = error("Cannot divide $AFF by $FAE")
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
# GenericNormExpr--UncAffExpr
(+){P,C,V<:Uncertain}(lhs::GenericNormExpr{P,C,V},rhs::UncAffExpr) = error("Cannot add $UNM to $FAE")
(-){P,C,V<:Uncertain}(lhs::GenericNormExpr{P,C,V},rhs::UncAffExpr) = error("Cannot subtract $UNM by $FAE")
(*){P,C,V<:Uncertain}(lhs::GenericNormExpr{P,C,V},rhs::UncAffExpr) = error("Cannot multiply $UNM by $FAE")
(/){P,C,V<:Uncertain}(lhs::GenericNormExpr{P,C,V},rhs::UncAffExpr) = error("Cannot divide $UNM by $FAE")


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
(-)(lhs::Uncertain, rhs::Variable) = UncAffExpr([rhs],[UncExpr(-1)],UncExpr(lhs))
(*)(lhs::Uncertain, rhs::Variable) = (*)(rhs, lhs)
(/)(lhs::Uncertain, rhs::Variable) = error("Cannot divide $UNC by a variable")
# Uncertain--AffExpr
(+)(lhs::Uncertain, rhs::AffExpr) = (+)(rhs, lhs)
(-)(lhs::Uncertain, rhs::AffExpr) = UncAffExpr(rhs.vars, map(UncExpr,-rhs.coeffs), UncExpr([lhs],[1.0],-rhs.constant))
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
# Uncertain--UncAffExpr
(+)(lhs::Uncertain, rhs::UncAffExpr) = UncAffExpr(rhs.vars,    rhs.coeffs,lhs+rhs.constant)
(-)(lhs::Uncertain, rhs::UncAffExpr) = UncAffExpr(rhs.vars,-1.*rhs.coeffs,lhs-rhs.constant)
(*)(lhs::Uncertain, rhs::UncAffExpr) = error("Cannot multiply $UNC by $FAE")
(/)(lhs::Uncertain, rhs::UncAffExpr) = error("Cannot divide $UNC by $FAE")

# UncExpr
# UncExpr--Number        - handled by JuMP
# UncExpr--Variable
(+)(lhs::UncExpr, rhs::Variable) = (+)(rhs,lhs)
(-)(lhs::UncExpr, rhs::Variable) = UncAffExpr([rhs],[UncExpr(-1)],lhs)
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
# UncExpr--UncAffExpr
(+)(lhs::UncExpr, rhs::UncAffExpr) = UncAffExpr(rhs.vars,    rhs.coeffs,lhs+rhs.constant)
(-)(lhs::UncExpr, rhs::UncAffExpr) = UncAffExpr(rhs.vars,-1.*rhs.coeffs,lhs-rhs.constant)
(*)(lhs::UncExpr, rhs::UncAffExpr) = error("Cannot multiply $UAE by $FAE")
(/)(lhs::UncExpr, rhs::UncAffExpr) = error("Cannot divide $UAE by $FAE")

# UncAffExpr
# UncAffExpr--Number     - handled by JuMP
# UncAffExpr--Variable
(+)(lhs::UncAffExpr, rhs::Variable) = (+)(rhs,lhs)
(-)(lhs::UncAffExpr, rhs::Variable) = UncAffExpr(vcat(lhs.vars,rhs),vcat(lhs.coeffs,UncExpr(-1)), lhs.constant)
(*)(lhs::UncAffExpr, rhs::Variable) = error("Cannot multiply $FAE by a variable")
(/)(lhs::UncAffExpr, rhs::Variable) = error("Cannot divide $FAE by a variable")
# UncAffExpr--AffExpr
(+)(lhs::UncAffExpr, rhs::AffExpr) = (+)(rhs,lhs)
(-)(lhs::UncAffExpr, rhs::AffExpr) = UncAffExpr(
  vcat(lhs.vars,    rhs.vars),
  vcat(lhs.coeffs, -rhs.coeffs),
  lhs.constant - rhs.constant)
(*)(lhs::UncAffExpr, rhs::AffExpr) = error("Cannot multiply $FAE by $AFF")
(/)(lhs::UncAffExpr, rhs::AffExpr) = error("Cannot divide $FAE by $AFF")
# UncAffExpr--GenericNormExpr
(+){P,C,V<:Uncertain}(lhs::UncAffExpr, rhs::GenericNormExpr{P,C,V}) = error("Cannot add $FAE to $UNM")
(-){P,C,V<:Uncertain}(lhs::UncAffExpr, rhs::GenericNormExpr{P,C,V}) = error("Cannot subtract $FAE by $UNM")
(*){P,C,V<:Uncertain}(lhs::UncAffExpr, rhs::GenericNormExpr{P,C,V}) = error("Cannot multiply $FAE by $UNM")
(/){P,C,V<:Uncertain}(lhs::UncAffExpr, rhs::GenericNormExpr{P,C,V}) = error("Cannot divide $FAE by $UNM")
# UncAffExpr--Uncertain
(+)(lhs::UncAffExpr, rhs::Uncertain) = UncAffExpr(lhs.vars,lhs.coeffs,lhs.constant+rhs)
(-)(lhs::UncAffExpr, rhs::Uncertain) = UncAffExpr(lhs.vars,lhs.coeffs,lhs.constant-rhs)
(*)(lhs::UncAffExpr, rhs::Uncertain) = error("Cannot multiply $FAE by $UNC")
(/)(lhs::UncAffExpr, rhs::Uncertain) = error("Cannot divide $FAE by $UNC")
# UncAffExpr--UncExpr
(+)(lhs::UncAffExpr, rhs::UncExpr) = UncAffExpr(lhs.vars,lhs.coeffs,lhs.constant+rhs)
(-)(lhs::UncAffExpr, rhs::UncExpr) = UncAffExpr(lhs.vars,lhs.coeffs,lhs.constant-rhs)
(*)(lhs::UncAffExpr, rhs::UncExpr) = error("Cannot multiply $FAE by $UAE")
(/)(lhs::UncAffExpr, rhs::UncExpr) = error("Cannot divide $FAE by $UAE")
# UncAffExpr--UncAffExpr
(*)(lhs::UncAffExpr, rhs::UncAffExpr) = error("Cannot multiply $FAE by $FAE")
(/)(lhs::UncAffExpr, rhs::UncAffExpr) = error("Cannot divide $FAE by $FAE")

#-----------------------------------------------------------------------
# For matrix operations
# 1. Number
# 2. Variable
# 3. [Generic]Norm
# 4. [Generic]AffExpr
# 5. QuadExpr <- We don't support any interactions with QuadExpr
# 6. [Generic]NormExpr
# 7. Uncertain
# 8. UncExpr
# 9. UncAffExpr
Base.promote_rule{R<:Real}(::Type{R},         ::Type{Uncertain} ) = UncExpr
Base.promote_rule{R<:Real}(::Type{R},         ::Type{UncExpr}   ) = UncExpr
Base.promote_rule{R<:Real}(::Type{R},         ::Type{UncAffExpr}) = UncAffExpr

Base.promote_rule{R<:Real}(::Type{Uncertain}, ::Type{R})          = UncExpr
Base.promote_rule(         ::Type{Uncertain}, ::Type{Uncertain})  = UncExpr
Base.promote_rule(         ::Type{Uncertain}, ::Type{UncExpr})    = UncExpr
Base.promote_rule(         ::Type{Uncertain}, ::Type{UncAffExpr}) = UncAffExpr

Base.promote_rule{R<:Real}(::Type{UncExpr},   ::Type{R})          = UncExpr
Base.promote_rule(         ::Type{UncExpr},   ::Type{Uncertain})  = UncExpr
Base.promote_rule(         ::Type{UncExpr},   ::Type{UncExpr})    = UncExpr
Base.promote_rule(         ::Type{UncExpr},   ::Type{UncAffExpr}) = UncAffExpr

Base.promote_rule{R<:Real}(::Type{UncAffExpr},::Type{R})          = UncAffExpr
Base.promote_rule(         ::Type{UncAffExpr},::Type{Uncertain})  = UncAffExpr
Base.promote_rule(         ::Type{UncAffExpr},::Type{UncExpr})    = UncAffExpr
Base.promote_rule(         ::Type{UncAffExpr},::Type{UncAffExpr}) = UncAffExpr


#-----------------------------------------------------------------------
# High-level operators like sum and dot are handled by JuMP
if VERSION < v"0.4-"
    import JuMP.vecdot
    export vecdot
end