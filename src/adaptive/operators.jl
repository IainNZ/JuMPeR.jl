#-----------------------------------------------------------------------
# JuMPeR  --  JuMP Extension for Robust Optimization
# http://github.com/IainNZ/JuMPeR.jl
#-----------------------------------------------------------------------
# Copyright (c) 2015: Iain Dunning
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
#-----------------------------------------------------------------------
# src/adaptive/macro.jl
# Adaptive robust optimization support - operators
# 1. Number
# 2. Variable
# 3. [Generic]Norm
# 4. [Generic]AffExpr   [Number * Variable]
# 5. QuadExpr <- We don't support any interactions with QuadExpr
# 6. [Generic]NormExpr
# 7. Uncertain
# 8. UncExpr            [Number * Uncertain]
# 9. UncAffExpr         [UncExpr * (Variable,Adaptive)]
# 10. Adaptive
# 11. AdaptAffExpr      [Number * Adaptive]
#-----------------------------------------------------------------------
    # Tab hint

#-----------------------------------------------------------------------
(+)(c::Number, x::Adaptive) = AdaptAffExpr(Adaptive[x], Float64[+1], c)
(-)(c::Number, x::Adaptive) = AdaptAffExpr(Adaptive[x], Float64[-1], c)
(*)(c::Number, x::Adaptive) = AdaptAffExpr(Adaptive[x], Float64[ c], 0)
#(/)(c::Number, x::Adaptive) = error()

#-----------------------------------------------------------------------

(+)(v::Variable, x::Adaptive) = UncAffExpr([v,x], [1,+1], 0)
(-)(v::Variable, x::Adaptive) = UncAffExpr([v,x], [1,-1], 0)
#(*)(v::Variable, x::Adaptive) = error()
#(/)(v::Variable, x::Adaptive) = error()

#-----------------------------------------------------------------------

(+)(a::AffExpr, x::Adaptive) = UncAffExpr(vcat(a.vars, x), vcat(a.coeffs, 1), a.constant)
(-)(a::AffExpr, x::Adaptive) = UncAffExpr(vcat(a.vars, x), vcat(a.coeffs,-1), a.constant)
#(*)(a::AffExpr, x::Adaptive) = error()
#(/)(a::AffExpr, x::Adaptive) = error()

(+)(a::AffExpr, b::AdaptAffExpr) = AdaptAffExpr(vcat(a.vars,   b.vars),
                                                vcat(a.coeffs, b.coeffs),
                                                a.constant + b.constant)
(-)(a::AffExpr, b::AdaptAffExpr) = AdaptAffExpr(vcat(a.vars,   b.vars),
                                                vcat(a.coeffs,-b.coeffs),
                                                a.constant - b.constant)
#(*)(a::AffExpr, b::AdaptAffExpr) = error()
#(/)(a::AffExpr, b::AdaptAffExpr) = error()

#-----------------------------------------------------------------------

(+)(u::Uncertain, x::Adaptive) = UncAffExpr([x], [ 1], u)
(-)(u::Uncertain, x::Adaptive) = UncAffExpr([x], [-1], u)
(*)(u::Uncertain, x::Adaptive) = UncAffExpr([x], [ u], 0)
#(/)(u::Uncertain, x::Adaptive) = error()

(+)(u::Uncertain, x::AdaptAffExpr) = UncAffExpr(copy(x.vars), copy(x.coeffs), u + x.constant)
(-)(u::Uncertain, x::AdaptAffExpr) = UncAffExpr(copy(x.vars),     -x.coeffs , u - x.constant)
(*)(u::Uncertain, x::AdaptAffExpr) = UncAffExpr(copy(x.vars),   u.*x.coeffs , u * x.constant)
#(/)(u::Uncertain, x::Adaptive) = error()

#-----------------------------------------------------------------------

(+)(u::UncExpr, x::Adaptive) = UncAffExpr([x], [ 1], u)
(-)(u::UncExpr, x::Adaptive) = UncAffExpr([x], [-1], u)
(*)(u::UncExpr, x::Adaptive) = UncAffExpr([x], [ u], 0)
#(/)(u::Uncertain, x::Adaptive) = error()

(+)(u::UncExpr, x::AdaptAffExpr) = UncAffExpr(copy(x.vars), copy(x.coeffs), u + x.constant)
(-)(u::UncExpr, x::AdaptAffExpr) = UncAffExpr(copy(x.vars),     -x.coeffs , u - x.constant)
(*)(u::UncExpr, x::AdaptAffExpr) = UncAffExpr(copy(x.vars),   u.*x.coeffs , u * x.constant)
#(/)(u::Uncertain, x::Adaptive) = error()

#-----------------------------------------------------------------------

+(a::UncAffExpr, b::AdaptAffExpr) = UncAffExpr(vcat(a.vars, b.vars),
                                               vcat(a.coeffs, map(UncExpr, b.coeffs)),
                                               a.constant + b.constant)
-(a::UncAffExpr, b::AdaptAffExpr) = UncAffExpr(vcat(a.vars, b.vars),
                                               vcat(a.coeffs, map(UncExpr,-b.coeffs)),
                                               a.constant - b.constant)
#*(a::UncAffExpr, b::AdaptAffExpr) = error()
#/(a::UncAffExpr, b::AdaptAffExpr) = error()

#-----------------------------------------------------------------------

(+)(x::Adaptive, c::Number) = AdaptAffExpr(Adaptive[x], Float64[  1],  c)
(-)(x::Adaptive, c::Number) = AdaptAffExpr(Adaptive[x], Float64[  1], -c)
(*)(x::Adaptive, c::Number) = AdaptAffExpr(Adaptive[x], Float64[  c],  0)
(/)(x::Adaptive, c::Number) = AdaptAffExpr(Adaptive[x], Float64[1/c],  0)

(+)(x::Adaptive, v::Variable) = UncAffExpr([x,v], [1,+1], 0)
(-)(x::Adaptive, v::Variable) = UncAffExpr([x,v], [1,-1], 0)
#(*)(x::Adaptive, v::Variable) = error()
#(/)(x::Adaptive, v::Variable) = error()

(+)(x::Adaptive, a::AffExpr) = UncAffExpr(vcat(x, a.vars), vcat(1, a.coeffs), a.constant)
(-)(x::Adaptive, a::AffExpr) = UncAffExpr(vcat(x, a.vars), vcat(1,-a.coeffs),-a.constant)
#(*)(x::Adaptive, a::AffExpr) = error()
#(/)(x::Adaptive, a::AffExpr) = error()

(+)(a::Adaptive, b::Adaptive) = AdaptAffExpr(Adaptive[a,b], Float64[1, 1], 0)
(-)(a::Adaptive, b::Adaptive) = AdaptAffExpr(Adaptive[a,b], Float64[1,-1], 0)
#(*)(a::Adaptive, b::Adaptive) = error()
#(/)(a::Adaptive, b::Adaptive) = error()

#-----------------------------------------------------------------------

(+)(a::AdaptAffExpr, b::UncExpr) = UncAffExpr(copy(a.vars),
                                              map(UncExpr,a.coeffs),
                                              copy(b))
