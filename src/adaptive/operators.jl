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
#   3. [Generic]Norm
# 4. [Generic]AffExpr   [Number * Variable]
#   5. QuadExpr <- We don't support any interactions with QuadExpr
#   6. [Generic]NormExpr
# 7. Uncertain
# 8. UncExpr            [Number * Uncertain]
# 9. UncAffExpr         [UncExpr * (Variable,Adaptive)]
# 10. Adaptive
# 11. AdaptAffExpr      [Number * Adaptive]
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

(+)(x::Adaptive, v::Variable) = AdaptAffExpr([x,v], [1,+1], 0)
(-)(x::Adaptive, v::Variable) = AdaptAffExpr([x,v], [1,-1], 0)
(*)(x::Adaptive, v::Variable) = error("Adaptive * Variable")
(/)(x::Adaptive, v::Variable) = error("Adaptive / Variable")

(+)(x::Adaptive, a::AffExpr) = AdaptAffExpr(vcat(x, a.vars), vcat(1,  a.coeffs),  a.constant)
(-)(x::Adaptive, a::AffExpr) = AdaptAffExpr(vcat(x, a.vars), vcat(1, -a.coeffs), -a.constant)
(*)(x::Adaptive, a::AffExpr) = error("Adaptive * AffExpr")
(/)(x::Adaptive, a::AffExpr) = error("Adaptive / AffExpr")

(+)(x::Adaptive, u::Uncertain) = UncAffExpr([x], [1],  u)
(-)(x::Adaptive, u::Uncertain) = UncAffExpr([x], [1], -u)
(*)(x::Adaptive, u::Uncertain) = UncAffExpr([x], [u],  0)
(/)(x::Adaptive, u::Uncertain) = error("Adaptive / Uncertain")

(+)(x::Adaptive, u::UncExpr) = UncAffExpr([x], [1],  u)
(-)(x::Adaptive, u::UncExpr) = UncAffExpr([x], [1], -u)
(*)(x::Adaptive, u::UncExpr) = UncAffExpr([x], [u],  0)
(/)(x::Adaptive, u::UncExpr) = error("Adaptive / UncExpr")

(+)(x::Adaptive, u::UncAffExpr) = UncAffExpr(vcat(x, u.vars), vcat(1,  u.coeffs),  u.constant)
(-)(x::Adaptive, u::UncAffExpr) = UncAffExpr(vcat(x, u.vars), vcat(1, -u.coeffs), -u.constant)
(*)(x::Adaptive, u::UncAffExpr) = error("Adaptive * UncAffExpr")
(/)(x::Adaptive, u::UncAffExpr) = error("Adaptive / UncAffExpr")

(+)(a::Adaptive, b::Adaptive) = AdaptAffExpr(Adaptive[a,b], Float64[1, 1], 0)
(-)(a::Adaptive, b::Adaptive) = AdaptAffExpr(Adaptive[a,b], Float64[1,-1], 0)
(*)(a::Adaptive, b::Adaptive) = error("Adaptive * Adaptive")
(/)(a::Adaptive, b::Adaptive) = error("Adaptive / Adaptive")

(+)(a::Adaptive, b::AdaptAffExpr) = AdaptAffExpr(vcat(a, b.vars), vcat(1,  b.coeffs),  b.constant)
(-)(a::Adaptive, b::AdaptAffExpr) = AdaptAffExpr(vcat(a, b.vars), vcat(1, -b.coeffs), -b.constant)
(*)(a::Adaptive, b::AdaptAffExpr) = error("Adaptive * AdaptAffExpr")
(/)(a::Adaptive, b::AdaptAffExpr) = error("Adaptive / AdaptAffExpr")

#-----------------------------------------------------------------------

(+)(a::AdaptAffExpr, c::Number) = AdaptAffExpr(copy(a.vars), copy(a.coeffs), a.constant + c)
(-)(a::AdaptAffExpr, c::Number) = AdaptAffExpr(copy(a.vars), copy(a.coeffs), a.constant - c)
(*)(a::AdaptAffExpr, c::Number) = AdaptAffExpr(copy(a.vars), copy(a.coeffs)*c, a.constant*c)
(/)(a::AdaptAffExpr, c::Number) = AdaptAffExpr(copy(a.vars), copy(a.coeffs)/c, a.constant/c)

(+)(x::AdaptAffExpr, v::Variable) = AdaptAffExpr(vcat(x.vars, v), vcat(x.coeffs,  1), x.constant)
(-)(x::AdaptAffExpr, v::Variable) = AdaptAffExpr(vcat(x.vars, v), vcat(x.coeffs, -1), x.constant)
(*)(x::AdaptAffExpr, v::Variable) = error("AdaptAffExpr * Variable")
(/)(x::AdaptAffExpr, v::Variable) = error("AdaptAffExpr / Variable")

(+)(x::AdaptAffExpr, a::AffExpr) = AdaptAffExpr(vcat(x.vars, a.vars), vcat(x.coeffs,  a.coeffs), x.constant + a.constant)
(-)(x::AdaptAffExpr, a::AffExpr) = AdaptAffExpr(vcat(x.vars, a.vars), vcat(x.coeffs, -a.coeffs), x.constant - a.constant)
(*)(x::AdaptAffExpr, a::AffExpr) = error("AdaptAffExpr * Variable")
(/)(x::AdaptAffExpr, a::AffExpr) = error("AdaptAffExpr / Variable")

(+)(a::AdaptAffExpr, u::Uncertain) = UncAffExpr(copy(a.vars), copy(a.coeffs), a.constant + u)
(-)(a::AdaptAffExpr, u::Uncertain) = UncAffExpr(copy(a.vars), copy(a.coeffs), a.constant - u)
(*)(a::AdaptAffExpr, u::Uncertain) = UncAffExpr(copy(a.vars), copy(a.coeffs) * u, a.constant * u)
(/)(a::AdaptAffExpr, u::Uncertain) = error("AdaptAffExpr / Uncertain")

(+)(a::AdaptAffExpr, b::UncExpr) = +( b, a)
(-)(a::AdaptAffExpr, b::UncExpr) = +(-b, a)
(*)(a::AdaptAffExpr, b::UncExpr) = *( b, a)
(/)(a::AdaptAffExpr, b::UncExpr) = /( b, a)  # not quite

(+)(a::AdaptAffExpr, b::UncAffExpr) = UncAffExpr(vcat(a.vars, b.vars),
                                                 vcat(map(UncExpr,a.coeffs), b.coeffs),
                                                 a.constant + b.constant)
(-)(a::AdaptAffExpr, b::UncAffExpr) = UncAffExpr(vcat(a.vars, b.vars),
                                                 vcat(map(UncExpr,a.coeffs), -b.coeffs),
                                                 a.constant - b.constant)
(*)(a::AdaptAffExpr, b::UncAffExpr) = error("AdaptAffExpr * UncAffExpr")
(/)(a::AdaptAffExpr, b::UncAffExpr) = error("AdaptAffExpr / UncAffExpr")

(+)(a::AdaptAffExpr, b::Adaptive) = AdaptAffExpr(vcat(a.vars,b), vcat(a.coeffs, 1), a.constant)
(-)(a::AdaptAffExpr, b::Adaptive) = AdaptAffExpr(vcat(a.vars,b), vcat(a.coeffs,-1), a.constant)
(*)(a::AdaptAffExpr, b::Adaptive) = error("AdaptAffExpr * Adaptive")
(/)(a::AdaptAffExpr, b::Adaptive) = error("AdaptAffExpr / Adaptive")

(+)(a::AdaptAffExpr, b::AdaptAffExpr) = AdaptAffExpr(vcat(a, b.vars), vcat(1,  b.coeffs), a.constant+b.constant)
(-)(a::AdaptAffExpr, b::AdaptAffExpr) = AdaptAffExpr(vcat(a, b.vars), vcat(1, -b.coeffs), a.constant-b.constant)
(*)(a::AdaptAffExpr, b::AdaptAffExpr) = error("AdaptAffExpr * AdaptAffExpr")
(/)(a::AdaptAffExpr, b::AdaptAffExpr) = error("AdaptAffExpr / AdaptAffExpr")
