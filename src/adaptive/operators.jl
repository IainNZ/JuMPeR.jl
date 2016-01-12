#-----------------------------------------------------------------------
# JuMPeR  --  JuMP Extension for Robust Optimization
# http://github.com/IainNZ/JuMPeR.jl
#-----------------------------------------------------------------------
# Copyright (c) 2016: Iain Dunning
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
    # Tab hint
#-----------------------------------------------------------------------

#-----------------------------------------------------------------------
(+)(c::Number, x::Adaptive) = AdaptAffExpr(Adaptive[x], Float64[+1], c)
(-)(c::Number, x::Adaptive) = AdaptAffExpr(Adaptive[x], Float64[-1], c)
(*)(c::Number, x::Adaptive) = AdaptAffExpr(Adaptive[x], Float64[ c], 0)
(/)(c::Number, x::Adaptive) = error("Number / Adaptive")

(+)(c::Number, x::AdaptAffExpr) = AdaptAffExpr(copy(x.vars), copy(x.coeffs), c + x.constant)
(-)(c::Number, x::AdaptAffExpr) = AdaptAffExpr(copy(x.vars),     -x.coeffs , c - x.constant)
(*)(c::Number, x::AdaptAffExpr) = AdaptAffExpr(copy(x.vars),  c * x.coeffs , c * x.constant)
(/)(c::Number, x::AdaptAffExpr) = error("Number / AdaptAffExpr")

#-----------------------------------------------------------------------

(+)(v::Variable, x::Adaptive) = AdaptAffExpr([v,x], [1,+1], 0)
(-)(v::Variable, x::Adaptive) = AdaptAffExpr([v,x], [1,-1], 0)
(*)(v::Variable, x::Adaptive) = error("Variable * Adaptive")
(/)(v::Variable, x::Adaptive) = error("Variable / Adaptive")

(+)(v::Variable, x::AdaptAffExpr) = AdaptAffExpr(vcat(v, x.vars), vcat(1,  x.coeffs),  x.constant)
(-)(v::Variable, x::AdaptAffExpr) = AdaptAffExpr(vcat(v, x.vars), vcat(1, -x.coeffs), -x.constant)
(*)(v::Variable, x::AdaptAffExpr) = error("Variable * AdaptAffExpr")
(/)(v::Variable, x::AdaptAffExpr) = error("Variable / AdaptAffExpr")

#-----------------------------------------------------------------------

(+)(a::AffExpr, x::Adaptive) = AdaptAffExpr(vcat(a.vars, x), vcat(a.coeffs,  1), a.constant)
(-)(a::AffExpr, x::Adaptive) = AdaptAffExpr(vcat(a.vars, x), vcat(a.coeffs, -1), a.constant)
(*)(a::AffExpr, x::Adaptive) = error("AffExpr * Adaptive")
(/)(a::AffExpr, x::Adaptive) = error("AffExpr / Adaptive")

(+)(a::AffExpr, b::AdaptAffExpr) = AdaptAffExpr(vcat(a.vars,   b.vars),
                                                vcat(a.coeffs, b.coeffs),
                                                a.constant + b.constant)
(-)(a::AffExpr, b::AdaptAffExpr) = AdaptAffExpr(vcat(a.vars,   b.vars),
                                                vcat(a.coeffs,-b.coeffs),
                                                a.constant - b.constant)
(*)(a::AffExpr, b::AdaptAffExpr) = error("AffExpr * AdaptAffExpr")
(/)(a::AffExpr, b::AdaptAffExpr) = error("AffExpr / AdaptAffExpr")

#-----------------------------------------------------------------------

(+)(u::Uncertain, x::Adaptive) = UncAffExpr([x], [ 1], u)
(-)(u::Uncertain, x::Adaptive) = UncAffExpr([x], [-1], u)
(*)(u::Uncertain, x::Adaptive) = UncAffExpr([x], [ u], 0)
(/)(u::Uncertain, x::Adaptive) = error("Uncertain / Adaptive")

(+)(u::Uncertain, x::AdaptAffExpr) = UncAffExpr(copy(x.vars), copy(x.coeffs), u + x.constant)
(-)(u::Uncertain, x::AdaptAffExpr) = UncAffExpr(copy(x.vars),     -x.coeffs , u - x.constant)
(*)(u::Uncertain, x::AdaptAffExpr) = UncAffExpr(copy(x.vars), u .* x.coeffs , u * x.constant)
(/)(u::Uncertain, x::AdaptAffExpr) = error("Uncertain / AdaptAffExpr")

#-----------------------------------------------------------------------

(+)(u::UncExpr, x::Adaptive) = UncAffExpr([x], [ 1], u)
(-)(u::UncExpr, x::Adaptive) = UncAffExpr([x], [-1], u)
(*)(u::UncExpr, x::Adaptive) = UncAffExpr([x], [ u], 0)
(/)(u::UncExpr, x::Adaptive) = error("UncExpr / Adaptive")

(+)(u::UncExpr, x::AdaptAffExpr) = UncAffExpr(copy(x.vars), copy(x.coeffs), u + x.constant)
(-)(u::UncExpr, x::AdaptAffExpr) = UncAffExpr(copy(x.vars),     -x.coeffs , u - x.constant)
(*)(u::UncExpr, x::AdaptAffExpr) = UncAffExpr(copy(x.vars), u .* x.coeffs , u * x.constant)
(/)(u::UncExpr, x::AdaptAffExpr) = error("UncExpr / AdaptAffExpr")

#-----------------------------------------------------------------------

+(a::UncAffExpr, x::Adaptive) = UncAffExpr(vcat(a.vars, x),
                                           vcat(a.coeffs, UncExpr(1)),
                                           a.constant)
-(a::UncAffExpr, x::Adaptive) = UncAffExpr(vcat(a.vars, x),
                                           vcat(a.coeffs, UncExpr(-1)),
                                           a.constant)
*(a::UncAffExpr, x::Adaptive) = error("UncAffExpr * Adaptive")
/(a::UncAffExpr, x::Adaptive) = error("UncAffExpr / Adaptive")

+(a::UncAffExpr, b::AdaptAffExpr) = UncAffExpr(vcat(a.vars, b.vars),
                                               vcat(a.coeffs, map(UncExpr, b.coeffs)),
                                               a.constant + b.constant)
-(a::UncAffExpr, b::AdaptAffExpr) = UncAffExpr(vcat(a.vars, b.vars),
                                               vcat(a.coeffs, map(UncExpr,-b.coeffs)),
                                               a.constant - b.constant)
*(a::UncAffExpr, b::AdaptAffExpr) = error("UncAffExpr * AdaptAffExpr")
/(a::UncAffExpr, b::AdaptAffExpr) = error("UncAffExpr / AdaptAffExpr")

#-----------------------------------------------------------------------

(+)(x::Adaptive, c::Number) = +(  c, x)
(-)(x::Adaptive, c::Number) = +( -c, x)
(*)(x::Adaptive, c::Number) = *(  c, x)
(/)(x::Adaptive, c::Number) = *(1/c, x)

(+)(x::Adaptive, v::Variable) = +( v, x)
(-)(x::Adaptive, v::Variable) = +(-v, x)
(*)(x::Adaptive, v::Variable) = *( v, x)
(/)(x::Adaptive, v::Variable) = /( v, x)  # not quite

(+)(x::Adaptive, a::AffExpr) = +( a, x)
(-)(x::Adaptive, a::AffExpr) = +(-a, x)
(*)(x::Adaptive, a::AffExpr) = *( a, x)
(/)(x::Adaptive, a::AffExpr) = /( a, x)  # not quite

(+)(x::Adaptive, u::Uncertain) = +( u, x)
(-)(x::Adaptive, u::Uncertain) = +(-u, x)
(*)(x::Adaptive, u::Uncertain) = *( u, x)
(/)(x::Adaptive, u::Uncertain) = /( u, x)  # not quite

(+)(a::Adaptive, b::Adaptive) = AdaptAffExpr(Adaptive[a,b], Float64[1, 1], 0)
(-)(a::Adaptive, b::Adaptive) = AdaptAffExpr(Adaptive[a,b], Float64[1,-1], 0)
(*)(a::Adaptive, b::Adaptive) = error("Adaptive * Adaptive")
(/)(a::Adaptive, b::Adaptive) = error("Adaptive / Adaptive")

#-----------------------------------------------------------------------

(+)(a::AdaptAffExpr, u::Uncertain) = +( u, a)
(-)(a::AdaptAffExpr, u::Uncertain) = +(-u, a)
(*)(a::AdaptAffExpr, u::Uncertain) = *( u, a)
(/)(a::AdaptAffExpr, u::Uncertain) = /( u, a)  # not quite

(+)(a::AdaptAffExpr, b::UncExpr) = +( b, a)
(-)(a::AdaptAffExpr, b::UncExpr) = +(-b, a)
(*)(a::AdaptAffExpr, b::UncExpr) = *( b, a)
(/)(a::AdaptAffExpr, b::UncExpr) = /( b, a)  # not quite
