#-----------------------------------------------------------------------
# JuMPeR  --  JuMP Extension for Robust Optimization
# http://github.com/IainNZ/JuMPeR.jl
#-----------------------------------------------------------------------
# Copyright (c) 2016: Iain Dunning
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
        @test_throws ErrorException 2 / ax
    end

    @testset "AdaptAffExpr" begin
        @test string(2 + aaff) == "3.4 ax + 5.1"
        @test string(2 - aaff) == "-3.4 ax - 1.1"
        @test string(2 * aaff) == "6.8 ax + 6.2"
        @test_throws ErrorException 2 / aaff
    end
end

@testset "Variable--" begin
    @testset "Adaptive" begin
        @test string(x + ax) == "x + ax"
        @test string(x - ax) == "x - ax"
        @test_throws ErrorException x * ax
        @test_throws ErrorException x / ax
    end

    @testset "AdaptAffExpr" begin
        @test string(x + aaff) == "x + 3.4 ax + 3.1"
        @test string(x - aaff) == "x - 3.4 ax - 3.1"
        @test_throws ErrorException x * aaff
        @test_throws ErrorException x / aaff
    end
end

@testset "AffExpr--" begin
    # aff  = 7.1 * x + 2.5
    @testset "Adaptive" begin
        @test string(aff + ax) == "7.1 x + ax + 2.5"
        @test string(aff - ax) == "7.1 x - ax + 2.5"
        @test_throws ErrorException aff * ax
        @test_throws ErrorException aff / ax
    end

    # aff  = 7.1 * x + 2.5
    # aaff = 3.4 * ax + 3.1
    @testset "AdaptAffExpr" begin
        @test string(aff + aaff) == "7.1 x + 3.4 ax + 5.6"
        @test string(aff - aaff) == "7.1 x - 3.4 ax - 0.6000000000000001"
        @test_throws ErrorException aff * aaff
        @test_throws ErrorException aff / aaff
    end
end

@testset "Uncertain--" begin
    # aff  = 7.1 * x + 2.5
    @testset "Adaptive" begin
        @test string(b + ax) ==  "ax + b"
        @test string(b - ax) == "-ax + b"
        @test string(b * ax) == "b ax"
        @test_throws ErrorException b / ax
    end

    # aaff = 3.4 * ax + 3.1
    @testset "AdaptAffExpr" begin
        @test string(b + aaff) ==  "3.4 ax + b + 3.1"
        @test string(b - aaff) == "-3.4 ax + b - 3.1"
        @test string(b * aaff) == "(3.4 b) ax + 3.1 b"
        @test_throws ErrorException b / aaff
    end
end

@testset "UncExpr--" begin
    # uex  = 2.3 * b + 5.5
    @testset "Adaptive" begin
        @test string(uex + ax) ==  "ax + 2.3 b + 5.5"
        @test string(uex - ax) == "-ax + 2.3 b + 5.5"
        @test string(uex * ax) == "(2.3 b + 5.5) ax"
        @test_throws ErrorException uex / ax
    end

    # uex  = 2.3 * b + 5.5
    # aaff = 3.4 * ax + 3.1
    @testset "AdaptAffExpr" begin
        @test string(uex + aaff) ==  "3.4 ax + 2.3 b + 8.6"
        @test string(uex - aaff) == "-3.4 ax + 2.3 b + 2.4"
        @test string(uex * aaff) == "(7.819999999999999 b + 18.7) ax + 7.13 b + 17.05"
        @test_throws ErrorException uex / aaff
    end
end

@testset "UncAffExpr--" begin
    # uaex = (5b+1)x + (2c + 3)
    @testset "Adaptive" begin
        @test string(uaex + ax) == "(5 b + 1) x + ax + 2 c + 3"
        @test string(uaex - ax) == "(5 b + 1) x - ax + 2 c + 3"
        @test_throws ErrorException uaex * ax
        @test_throws ErrorException uaex / ax
    end

    # uaex = (5b+1)x + (2c + 3)
    # aaff = 3.4 * ax + 3.1
    @testset "AdaptAffExpr" begin
        @test string(uaex + aaff) == "(5 b + 1) x + 3.4 ax + 2 c + 6.1"
        @test string(uaex - aaff) == "(5 b + 1) x - 3.4 ax + 2 c - 0.10000000000000009"
        @test_throws ErrorException uaex * aaff
        @test_throws ErrorException uaex / aaff
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
        @test string(ax + y) == "y + ax"
        @test string(ax - y) == "-y + ax"
        @test_throws ErrorException ax * y
        @test_throws ErrorException ax / y
    end

    @testset "AffExpr" begin
        @test string(ax + aff) == "7.1 x + ax + 2.5"
        @test string(ax - aff) == "-7.1 x + ax - 2.5"
        @test_throws ErrorException ax * aff
        @test_throws ErrorException ax / aff
    end

    @testset "Uncertain" begin
        @test string(ax + b) == "ax + b"
        @test string(ax - b) == "ax + -b"
        @test string(ax * b) == "b ax"
        @test_throws ErrorException ax / b
    end

    @testset "UncExpr" begin
        @test string(ax + uex) == "ax + 2.3 b + 5.5"
        @test string(ax - uex) == "ax + -2.3 b - 5.5"
        @test string(ax * uex) == "(2.3 b + 5.5) ax"
        @test_throws ErrorException ax / uex
    end

    @testset "UncAffExpr" begin
        @test string(ax + uaex) == "ax + (5 b + 1) x + 2 c + 3"
        @test string(ax - uaex) == "ax + (-5 b - 1) x + -2 c - 3"
        @test_throws ErrorException ax * uaex
        @test_throws ErrorException ax / uaex
    end

    @testset "Adaptive" begin
        @test string(ax + ay) == "ax + ay"
        @test string(ax - ay) == "ax - ay"
        @test_throws ErrorException ax * ay
        @test_throws ErrorException ax / ay
    end

    # aaff = 3.4 * ax + 3.1
    @testset "AdaptAffExpr" begin
        @test string(ay + aaff) == "ay + 3.4 ax + 3.1"
        @test string(ay - aaff) == "ay - 3.4 ax - 3.1"
        @test_throws ErrorException ay * aaff
        @test_throws ErrorException ay / aaff
    end
end  # "Adaptive--"

@testset "AdaptAffExpr--" begin
    @testset "Number" begin
        @test string(aaff + 2) == "3.4 ax + 5.1"
        @test string(aaff - 2) == "3.4 ax + 1.1"
        @test string(aaff * 2) == "6.8 ax + 6.2"
        @test string(aaff / 2) == "1.7 ax + 1.55"
    end

    @testset "Variable" begin
        @test string(aaff + x) == "x + 3.4 ax + 3.1"
        @test string(aaff - x) == "-x + 3.4 ax + 3.1"
        @test_throws ErrorException aaff * x
        @test_throws ErrorException aaff / x
    end

    # aff  = 7.1 * x + 2.5
    # aaff = 3.4 * ax + 3.1
    @testset "AffExpr" begin
        @test string(aaff + aff) == "7.1 x + 3.4 ax + 5.6"
        @test string(aaff - aff) == "-7.1 x + 3.4 ax + 0.6000000000000001"
        @test_throws ErrorException aaff * aff
        @test_throws ErrorException aaff / aff
    end

    @testset "Uncertain" begin
        @test string(aaff + b) == "3.4 ax + b + 3.1"
        @test string(aaff - b) == "3.4 ax + -b + 3.1"
        @test string(aaff * b) == "(3.4 b) ax + 3.1 b"
        @test_throws ErrorException aaff / b
    end

    # aaff = 3.4 * ax + 3.1
    # uex  = 2.3 * b + 5.5
    @testset "UncExpr" begin
        @test string(aaff + uex) == "3.4 ax + 2.3 b + 8.6"
        @test string(aaff - uex) == "3.4 ax + -2.3 b - 2.4"
        @test string(aaff * uex) == "(7.819999999999999 b + 18.7) ax + 7.13 b + 17.05"
        @test_throws ErrorException aaff / uex
    end

    # aaff = 3.4 * ax + 3.1
    # uaex = (5b+1)x + (2c + 3)
    @testset "UncAffExpr" begin
        @test string(aaff + uaex) == "3.4 ax + (5 b + 1) x + 2 c + 6.1"
        @test string(aaff - uaex) == "3.4 ax + (-5 b - 1) x + -2 c + 0.10000000000000009"
        @test_throws ErrorException aaff * uaex
        @test_throws ErrorException aaff / uaex
    end

    @testset "Adaptive" begin
        @test string(aaff + ay) == "3.4 ax + ay + 3.1"
        @test string(aaff - ay) == "3.4 ax - ay + 3.1"
        @test_throws ErrorException aaff * ay
        @test_throws ErrorException aaff / ay
    end
end  # "AdaptAffExpr--"

end  # "adp_ops"
