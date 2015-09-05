#-----------------------------------------------------------------------
# JuMPeR  --  JuMP Extension for Robust Optimization
# http://github.com/IainNZ/JuMPeR.jl
#-----------------------------------------------------------------------
# Copyright (c) 2015: Iain Dunning
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
#-----------------------------------------------------------------------
# test/adp_ops.jl
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

using JuMP, JuMPeR
using BaseTestNext
import JuMP: REPLMode, IJuliaMode

m = RobustModel()
@defVar(m, x)
@defVar(m, y)
@defAdaptVar(m, ax)
@defAdaptVar(m, ay)
@defUnc(m, b)
@defUnc(m, c)

aff  = 7.1 * x + 2.5
aaff = 3.4 * ax + 3.1
uex  = 2.3 * b + 5.5
uaex = (5b+1)x + (2c + 3)

@testset "adp_ops" begin

@testset "Number--" begin
    @testset "Adaptive" begin
        @test string(2 + ax) == "ax + 2"
        @test string(2 - ax) == "-ax + 2"
        @test string(2 * ax) == "2 ax"
        @test_throws MethodError 2 / ax
    end

    @testset "AdaptAffExpr" begin
        @test string(2 + aaff) == "3.4 ax + 5.1"
        @test string(2 - aaff) == "-3.4 ax - 1.1"
        @test string(2 * aaff) == "6.8 ax + 6.2"
        @test_throws MethodError 2 / aaff
    end
end

@testset "Variable--" begin
    @testset "Adaptive" begin
        @test string(x + ax) == "x + ax"
        @test string(x - ax) == "x - ax"
        @test_throws MethodError x * ax
        @test_throws MethodError x / ax
    end

    @testset "AdaptAffExpr" begin
        @test string(x + aaff) == "x + 3.4 ax + 3.1"
        @test string(x - aaff) == "x - 3.4 ax - 3.1"
        @test_throws MethodError x * aaff
        @test_throws MethodError x / aaff
    end
end

@testset "AffExpr--" begin
    # aff  = 7.1 * x + 2.5
    @testset "Adaptive" begin
        @test string(aff + ax) == "7.1 x + ax + 2.5"
        @test string(aff - ax) == "7.1 x - ax + 2.5"
        @test_throws MethodError aff * ax
        @test_throws MethodError aff / ax
    end

    # aff  = 7.1 * x + 2.5
    # aaff = 3.4 * ax + 3.1
    @testset "AdaptAffExpr" begin
        @test string(aff + aaff) == "7.1 x + 3.4 ax + 5.6"
        @test string(aff - aaff) == "7.1 x - 3.4 ax - 0.6000000000000001"
        @test_throws MethodError aff * aaff
        @test_throws MethodError aff / aaff
    end
end

@testset "Uncertain--" begin
    # aff  = 7.1 * x + 2.5
    @testset "Adaptive" begin
        @test string(b + ax) ==  "ax + b"
        @test string(b - ax) == "-ax + b"
        @test string(b * ax) == "b ax"
        @test_throws MethodError b / ax
    end

    # aaff = 3.4 * ax + 3.1
    @testset "AdaptAffExpr" begin
        @test string(b + aaff) ==  "3.4 ax + b + 3.1"
        @test string(b - aaff) == "-3.4 ax + b - 3.1"
        @test string(b * aaff) == "(3.4 b) ax + 3.1 b"
        @test_throws MethodError b / aaff
    end
end

@testset "UncExpr--" begin
    # uex  = 2.3 * b + 5.5
    @testset "Adaptive" begin
        @test string(uex + ax) ==  "ax + 2.3 b + 5.5"
        @test string(uex - ax) == "-ax + 2.3 b + 5.5"
        @test string(uex * ax) == "(2.3 b + 5.5) ax"
        @test_throws MethodError uex / ax
    end

    # uex  = 2.3 * b + 5.5
    # aaff = 3.4 * ax + 3.1
    @testset "AdaptAffExpr" begin
        @test string(uex + aaff) ==  "3.4 ax + 2.3 b + 8.6"
        @test string(uex - aaff) == "-3.4 ax + 2.3 b + 2.4"
        @test string(uex * aaff) == "(7.819999999999999 b + 18.7) ax + 7.13 b + 17.05"
        @test_throws MethodError uex / aaff
    end
end

@testset "UncAffExpr--" begin
    # uaex = (5b+1)x + (2c + 3)
    @testset "Adaptive" begin
        @test string(uaex + ax) == "(5 b + 1) x + ax + 2 c + 3"
        @test string(uaex - ax) == "(5 b + 1) x - ax + 2 c + 3"
        @test_throws MethodError uaex * ax
        @test_throws MethodError uaex / ax
    end

    # uaex = (5b+1)x + (2c + 3)
    # aaff = 3.4 * ax + 3.1
    @testset "AdaptAffExpr" begin
        @test string(uaex + aaff) == "(5 b + 1) x + 3.4 ax + 2 c + 6.1"
        @test string(uaex - aaff) == "(5 b + 1) x - 3.4 ax + 2 c - 0.10000000000000009"
        @test_throws MethodError uaex * aaff
        @test_throws MethodError uaex / aaff
    end
end

@testset "Adaptive--" begin
    @testset "Number" begin
        @test string(ax + 2) == "ax + 2"
        @test string(ax - 2) == "ax - 2"
        @test string(ax * 2) == "2 ax"
        @test string(ax / 2) == "0.5 ax"
    end

    @testset "Variable" begin
        @test string(ax + y) == "ax + y"
        @test string(ax - y) == "ax - y"
        @test_throws MethodError ax * y
        @test_throws MethodError ax / y
    end

    @testset "AffExpr" begin
        @test string(ax + aff) == "ax + 7.1 x + 2.5"
        @test string(ax - aff) == "ax - 7.1 x + -2.5"
        @test_throws MethodError ax * aff
        @test_throws MethodError ax / aff
    end

    @testset "Adaptive" begin
        @test string(ax + ay) == "ax + ay"
        @test string(ax - ay) == "ax - ay"
        @test_throws MethodError ax * ay
        @test_throws MethodError ax / ay
    end
end

end
