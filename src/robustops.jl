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
# 8. UAffExpr
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
(+)(lhs::Number, rhs::Uncertain) = UAffExpr(  1, rhs, lhs)
(-)(lhs::Number, rhs::Uncertain) = UAffExpr( -1, rhs, lhs)
(*)(lhs::Number, rhs::Uncertain) = UAffExpr(lhs, rhs)
(/)(lhs::Number, rhs::Uncertain) = error("Cannot divide a number by $UNC")
# Number--UAffExpr        - handled by JuMP
# Number--FullAffExpr     - handled by JuMP

#-----------------------------------------------------------------------
# 2. Variable
# Variable--Uncertain
(+)(lhs::Variable, rhs::Uncertain) = FullAffExpr([lhs],[UAffExpr(  1)],UAffExpr( 1,rhs))
(-)(lhs::Variable, rhs::Uncertain) = FullAffExpr([lhs],[UAffExpr(  1)],UAffExpr(-1,rhs))
(*)(lhs::Variable, rhs::Uncertain) = FullAffExpr([lhs],[UAffExpr(rhs)],UAffExpr())
(/)(lhs::Variable, rhs::Uncertain) = error("Cannot divide a variable by $UNC")
# Variable--UAffExpr
(+)(lhs::Variable, rhs::UAffExpr) = FullAffExpr([lhs],[UAffExpr(1)],       rhs)
(-)(lhs::Variable, rhs::UAffExpr) = FullAffExpr([lhs],[UAffExpr(1)],      -rhs)
(*)(lhs::Variable, rhs::UAffExpr) = FullAffExpr([lhs],[        rhs],UAffExpr())
(/)(lhs::Variable, rhs::UAffExpr) = error("Cannot divide a variable by $UAE")
# Variable--FullAffExpr
(+)(lhs::Variable, rhs::FullAffExpr) = FullAffExpr(vcat(rhs.vars,lhs),vcat(    rhs.coeffs,UAffExpr(1)), rhs.constant)
(-)(lhs::Variable, rhs::FullAffExpr) = FullAffExpr(vcat(rhs.vars,lhs),vcat(-1.*rhs.coeffs,UAffExpr(1)),-rhs.constant)
(*)(lhs::Variable, rhs::FullAffExpr) = error("Cannot multiply a variable by $FAE")
(/)(lhs::Variable, rhs::FullAffExpr) = error("Cannot divide a variable by $FAE")

#-----------------------------------------------------------------------
# 3. Norm


#-----------------------------------------------------------------------
# 4. AffExpr
# AffExpr--Uncertain
(+)(lhs::AffExpr, rhs::Uncertain) = FullAffExpr(lhs.vars, map(UAffExpr,lhs.coeffs), UAffExpr( 1,rhs,lhs.constant))
(-)(lhs::AffExpr, rhs::Uncertain) = FullAffExpr(lhs.vars, map(UAffExpr,lhs.coeffs), UAffExpr(-1,rhs,lhs.constant))
(*)(lhs::AffExpr, rhs::Uncertain) = FullAffExpr(lhs.vars,         rhs.*lhs.coeffs , UAffExpr(lhs.constant,rhs))
(/)(lhs::AffExpr, rhs::Uncertain) = error("Cannot divide $AFF by $UNC")
# AffExpr-UAffExpr
(+)(lhs::AffExpr, rhs::UAffExpr)  = FullAffExpr(lhs.vars, map(UAffExpr,lhs.coeffs), lhs.constant+rhs)
(-)(lhs::AffExpr, rhs::UAffExpr)  = FullAffExpr(lhs.vars, map(UAffExpr,lhs.coeffs), lhs.constant-rhs)
(*)(lhs::AffExpr, rhs::UAffExpr)  = FullAffExpr(lhs.vars,         rhs.*lhs.coeffs , lhs.constant*rhs)
(/)(lhs::AffExpr, rhs::UAffExpr)  = error("Cannot divide $AFF by $UAE")
# AffExpr-FullAffExpr
(+)(lhs::AffExpr, rhs::FullAffExpr) = FullAffExpr(
  vcat(lhs.vars, rhs.vars),
  vcat(map(UAffExpr,lhs.coeffs), rhs.coeffs),
  lhs.constant + rhs.constant)
(-)(lhs::AffExpr, rhs::FullAffExpr) = FullAffExpr(
  vcat(lhs.vars, rhs.vars),
  vcat(map(UAffExpr,lhs.coeffs), -1.*rhs.coeffs),
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
# GenericNormExpr--UAffExpr
(+){P,C,V<:Uncertain}(lhs::GenericNormExpr{P,C,V},rhs::UAffExpr) = error("Cannot add $UNM to $UAE")
(-){P,C,V<:Uncertain}(lhs::GenericNormExpr{P,C,V},rhs::UAffExpr) = error("Cannot subtract $UNM by $UAE")
(*){P,C,V<:Uncertain}(lhs::GenericNormExpr{P,C,V},rhs::UAffExpr) = error("Cannot multiply $UNM by $UAE")
(/){P,C,V<:Uncertain}(lhs::GenericNormExpr{P,C,V},rhs::UAffExpr) = error("Cannot divide $UNM by $UAE")
# GenericNormExpr--FullAffExpr
(+){P,C,V<:Uncertain}(lhs::GenericNormExpr{P,C,V},rhs::FullAffExpr) = error("Cannot add $UNM to $FAE")
(-){P,C,V<:Uncertain}(lhs::GenericNormExpr{P,C,V},rhs::FullAffExpr) = error("Cannot subtract $UNM by $FAE")
(*){P,C,V<:Uncertain}(lhs::GenericNormExpr{P,C,V},rhs::FullAffExpr) = error("Cannot multiply $UNM by $FAE")
(/){P,C,V<:Uncertain}(lhs::GenericNormExpr{P,C,V},rhs::FullAffExpr) = error("Cannot divide $UNM by $FAE")


#-----------------------------------------------------------------------
# Uncertain
(-)(lhs::Uncertain) = UAffExpr(-1,lhs)
# Uncertain--Number
(+)(lhs::Uncertain, rhs::Number) = UAffExpr(    1, lhs,  rhs)
(-)(lhs::Uncertain, rhs::Number) = UAffExpr(    1, lhs, -rhs)
(*)(lhs::Uncertain, rhs::Number) = UAffExpr(  rhs, lhs)
(/)(lhs::Uncertain, rhs::Number) = UAffExpr(1/rhs, lhs)
# Uncertain--Variable
(+)(lhs::Uncertain, rhs::Variable) = (+)(rhs, lhs)
(-)(lhs::Uncertain, rhs::Variable) = FullAffExpr([rhs],[UAffExpr(-1)],UAffExpr(lhs))
(*)(lhs::Uncertain, rhs::Variable) = (*)(rhs, lhs)
(/)(lhs::Uncertain, rhs::Variable) = error("Cannot divide $UNC by a variable")
# Uncertain--AffExpr
(+)(lhs::Uncertain, rhs::AffExpr) = (+)(rhs, lhs)
(-)(lhs::Uncertain, rhs::AffExpr) = FullAffExpr(rhs.vars, map(UAffExpr,-rhs.coeffs), UAffExpr(1,lhs,-rhs.constant))
(*)(lhs::Uncertain, rhs::AffExpr) = (*)(rhs, lhs)
(/)(lhs::Uncertain, rhs::AffExpr) = error("Cannot divide $UNC by $AFF")
# Uncertain--GenericNormExpr
(+){P,C,V<:Uncertain}(lhs::Uncertain, rhs::GenericNormExpr{P,C,V}) = error("Cannot add $UNC to $UNM")
(-){P,C,V<:Uncertain}(lhs::Uncertain, rhs::GenericNormExpr{P,C,V}) = error("Cannot subtract $UNC by $UNM")
(*){P,C,V<:Uncertain}(lhs::Uncertain, rhs::GenericNormExpr{P,C,V}) = error("Cannot multiply $UNC by $UNM")
(/){P,C,V<:Uncertain}(lhs::Uncertain, rhs::GenericNormExpr{P,C,V}) = error("Cannot divide $UNC by $UNM")
# Uncertain--Uncertain
(+)(lhs::Uncertain, rhs::Uncertain) = UAffExpr([lhs,rhs], Float64[1, 1], 0.0)
(-)(lhs::Uncertain, rhs::Uncertain) = UAffExpr([lhs,rhs], Float64[1,-1], 0.0)
(*)(lhs::Uncertain, rhs::Uncertain) = error("Cannot multiply $UNC by $UNC")
(/)(lhs::Uncertain, rhs::Uncertain) = error("Cannot divide $UNC by $UNC")
# Uncertain--UAffExpr
(+)(lhs::Uncertain, rhs::UAffExpr) = UAffExpr(vcat(lhs,rhs.vars),vcat(1.0, rhs.coeffs), rhs.constant)
(-)(lhs::Uncertain, rhs::UAffExpr) = UAffExpr(vcat(lhs,rhs.vars),vcat(1.0,-rhs.coeffs),-rhs.constant)
(*)(lhs::Uncertain, rhs::UAffExpr) = error("Cannot multiply $UNC by $UAE")
(/)(lhs::Uncertain, rhs::UAffExpr) = error("Cannot divide $UNC by $UAE")
# Uncertain--FullAffExpr
(+)(lhs::Uncertain, rhs::FullAffExpr) = FullAffExpr(rhs.vars,    rhs.coeffs,lhs+rhs.constant)
(-)(lhs::Uncertain, rhs::FullAffExpr) = FullAffExpr(rhs.vars,-1.*rhs.coeffs,lhs-rhs.constant)
(*)(lhs::Uncertain, rhs::FullAffExpr) = error("Cannot multiply $UNC by $FAE")
(/)(lhs::Uncertain, rhs::FullAffExpr) = error("Cannot divide $UNC by $FAE")

# UAffExpr
# UAffExpr--Number        - handled by JuMP
# UAffExpr--Variable
(+)(lhs::UAffExpr, rhs::Variable) = (+)(rhs,lhs)
(-)(lhs::UAffExpr, rhs::Variable) = FullAffExpr([rhs],[UAffExpr(-1)],lhs)
(*)(lhs::UAffExpr, rhs::Variable) = (*)(rhs,lhs)
(/)(lhs::UAffExpr, rhs::Variable) = error("Cannot divide $UAE by a variable")
# UAffExpr--AffExpr
(+)(lhs::UAffExpr, rhs::AffExpr) = (+)( rhs,lhs)
(-)(lhs::UAffExpr, rhs::AffExpr) = (+)(-rhs,lhs)
(*)(lhs::UAffExpr, rhs::AffExpr) = (*)( rhs,lhs)
(/)(lhs::UAffExpr, rhs::AffExpr) = error("Cannot divide $UAE by $AFF")
# UAffExpr--GenericNormExpr
(+){P,C,V<:Uncertain}(lhs::UAffExpr, rhs::GenericNormExpr{P,C,V}) = error("Cannot add $UAE to $UNM")
(-){P,C,V<:Uncertain}(lhs::UAffExpr, rhs::GenericNormExpr{P,C,V}) = error("Cannot subtract $UAE by $UNM")
(*){P,C,V<:Uncertain}(lhs::UAffExpr, rhs::GenericNormExpr{P,C,V}) = error("Cannot multiply $UAE by $UNM")
(/){P,C,V<:Uncertain}(lhs::UAffExpr, rhs::GenericNormExpr{P,C,V}) = error("Cannot divide $UAE by $UNM")
# UAffExpr--Uncertain
(+)(lhs::UAffExpr, rhs::Uncertain) = (+)(rhs,lhs)
(-)(lhs::UAffExpr, rhs::Uncertain) = UAffExpr(vcat(rhs,lhs.vars),vcat(-1.0,lhs.coeffs),lhs.constant)
(*)(lhs::UAffExpr, rhs::Uncertain) = (*)(rhs,lhs)
(/)(lhs::UAffExpr, rhs::Uncertain) = error("Cannot divide $UAE by $UNC")
# UAffExpr--UAffExpr
(*)(lhs::UAffExpr, rhs::UAffExpr) = error("Cannot multiply $UAE by $UAE")
(/)(lhs::UAffExpr, rhs::UAffExpr) = error("Cannot divide $UAE by $UAE")
# UAffExpr--FullAffExpr
(+)(lhs::UAffExpr, rhs::FullAffExpr) = FullAffExpr(rhs.vars,    rhs.coeffs,lhs+rhs.constant)
(-)(lhs::UAffExpr, rhs::FullAffExpr) = FullAffExpr(rhs.vars,-1.*rhs.coeffs,lhs-rhs.constant)
(*)(lhs::UAffExpr, rhs::FullAffExpr) = error("Cannot multiply $UAE by $FAE")
(/)(lhs::UAffExpr, rhs::FullAffExpr) = error("Cannot divide $UAE by $FAE")

# FullAffExpr
# FullAffExpr--Number     - handled by JuMP
# FullAffExpr--Variable
(+)(lhs::FullAffExpr, rhs::Variable) = (+)(rhs,lhs)
(-)(lhs::FullAffExpr, rhs::Variable) = FullAffExpr(vcat(lhs.vars,rhs),vcat(lhs.coeffs,UAffExpr(-1)), lhs.constant)
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
# UAffExpr--GenericNormExpr
(+){P,C,V<:Uncertain}(lhs::FullAffExpr, rhs::GenericNormExpr{P,C,V}) = error("Cannot add $FAE to $UNM")
(-){P,C,V<:Uncertain}(lhs::FullAffExpr, rhs::GenericNormExpr{P,C,V}) = error("Cannot subtract $FAE by $UNM")
(*){P,C,V<:Uncertain}(lhs::FullAffExpr, rhs::GenericNormExpr{P,C,V}) = error("Cannot multiply $FAE by $UNM")
(/){P,C,V<:Uncertain}(lhs::FullAffExpr, rhs::GenericNormExpr{P,C,V}) = error("Cannot divide $FAE by $UNM")
# FullAffExpr--Uncertain
(+)(lhs::FullAffExpr, rhs::Uncertain) = FullAffExpr(lhs.vars,lhs.coeffs,lhs.constant+rhs)
(-)(lhs::FullAffExpr, rhs::Uncertain) = FullAffExpr(lhs.vars,lhs.coeffs,lhs.constant-rhs)
(*)(lhs::FullAffExpr, rhs::Uncertain) = error("Cannot multiply $FAE by $UNC")
(/)(lhs::FullAffExpr, rhs::Uncertain) = error("Cannot divide $FAE by $UNC")
# FullAffExpr--UAffExpr
(+)(lhs::FullAffExpr, rhs::UAffExpr) = FullAffExpr(lhs.vars,lhs.coeffs,lhs.constant+rhs)
(-)(lhs::FullAffExpr, rhs::UAffExpr) = FullAffExpr(lhs.vars,lhs.coeffs,lhs.constant-rhs)
(*)(lhs::FullAffExpr, rhs::UAffExpr) = error("Cannot multiply $FAE by $UAE")
(/)(lhs::FullAffExpr, rhs::UAffExpr) = error("Cannot divide $FAE by $UAE")
# FullAffExpr--FullAffExpr
(*)(lhs::FullAffExpr, rhs::FullAffExpr) = error("Cannot multiply $FAE by $FAE")
(/)(lhs::FullAffExpr, rhs::FullAffExpr) = error("Cannot divide $FAE by $FAE")


#-----------------------------------------------------------------------
# High-level operators like sum and dot are handled by JuMP
if VERSION < v"0.4-"
    import JuMP.vecdot
    export vecdot
end


#-----------------------------------------------------------------------
# Matrix operations
import JuMP: _multiply_type
# ----
_multiply_type{R<:Real}(::Type{R},::Type{Uncertain})    = UAffExpr
_multiply_type{R<:Real}(::Type{R},::Type{UAffExpr})     = UAffExpr
_multiply_type{R<:Real}(::Type{R},::Type{FullAffExpr})  = FullAffExprr
# ----
_multiply_type(::Type{Variable},::Type{Uncertain})      = FullAffExpr
_multiply_type(::Type{Variable},::Type{UAffExpr})       = FullAffExpr
# ----
_multiply_type(::Type{AffExpr}, ::Type{Uncertain})      = UAffExpr
_multiply_type(::Type{AffExpr}, ::Type{UAffExpr})       = FullAffExpr
_multiply_type(::Type{AffExpr}, ::Type{FullAffExpr})    = FullAffExpr
# ----
_multiply_type{R<:Real}(::Type{Uncertain},::Type{R})    = UAffExpr
_multiply_type(::Type{Uncertain},::Type{Variable})      = FullAffExpr
_multiply_type(::Type{Uncertain},::Type{AffExpr})       = FullAffExpr
# ----
_multiply_type{R<:Real}(::Type{UAffExpr},::Type{R})     = UAffExpr
_multiply_type(::Type{UAffExpr},::Type{Variable})       = FullAffExpr
_multiply_type(::Type{UAffExpr},::Type{AffExpr})        = FullAffExpr
# ----
_multiply_type{R<:Real}(::Type{FullAffExpr},::Type{R})  = FullAffExpr