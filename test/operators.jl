#-----------------------------------------------------------------------
# JuMPeR  --  JuMP Extension for Robust Optimization
# http://github.com/IainNZ/JuMPeR.jl
#-----------------------------------------------------------------------
# Copyright (c) 2016: Iain Dunning
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
#-----------------------------------------------------------------------
# test/operators.jl
# Testing for all operator overloads and norm construction
#-----------------------------------------------------------------------

using JuMP, JuMPeR
using Test
import JuMP: REPLMode, IJuliaMode

# To ensure the tests work on both Windows and Linux/OSX, we need
# to use the correct comparison operators in the strings
const leq = JuMP.repl[:leq]
const geq = JuMP.repl[:geq]
const  eq = JuMP.repl[:eq]

@testset "Operators" begin
printstyled("Operators...\n", color = :yellow)

@testset "Basics" begin

# Setup
m = RobustModel()
@variable(m, w)
@variable(m, x)
@variable(m, y)
@variable(m, z)
@adaptive(m, ax)
@adaptive(m, ay)
@uncertain(m, a)
@uncertain(m, b)
@uncertain(m, c)

aff   = 7.1 * x + 2.5
aff2  = 1.2 * y + 1.2
uaff  = 2.3 * a + 5.5
uaff2 = 3.4 * b + 1.1
faff  = (5a+1)x + (2b + 3)
aaff = 3.4 * ax + 3.1
uex  = 2.3 * b + 5.5
uaex = (5b+1)x + (2c + 3)

@test sprint(print, a)  == "a"
@test sprint( show, a)  == "a"
@test string(aff)     == "7.1 x + 2.5"
@test string(aff2)    == "1.2 y + 1.2"
@test string(uaff)    == "2.3 a + 5.5"
@test string(uaff2)   == "3.4 b + 1.1"
@test string(faff)    == "(5 a + 1) x + 2 b + 3"
@test string(x + JuMPeR.UncExpr()) == sprint(print, x)


@testset "Number--... tests" begin
    # Number--Uncertain
    @test string(4.13 + a) == "a + 4.13"
    @test string(3.16 - a) == "-a + 3.16"
    @test string(5.23 * a) == "5.23 a"
    @test_throws ErrorException 2.94 / a
    # Number--UncExpr
    @test string(2.3 + uaff) == "2.3 a + 7.8"
    @test string(1.5 - uaff) == "-2.3 a - 4"
    @test string(2.0 * uaff) == "4.6 a + 11"
    @test_throws MethodError 2.94 / uaff
    # Number--UncVarExpr
    @test string(2.3 + faff) == "(5 a + 1) x + 2 b + 5.3"
    @test string(1.0 - faff) == "(-5 a - 1) x + -2 b - 2"
    @test string(2.0 * faff) == "(10 a + 2) x + 4 b + 6"
    @test_throws MethodError 2.94 / faff
    @testset "Adaptive" begin
        @test string(2 + ax) == "ax + 2"
        @test string(2 - ax) == "-ax + 2"
        @test string(2 * ax) == "2 ax"
        @test_throws ErrorException 2 / ax
    end
    @testset "AdaptExpr" begin
        @test string(2 + aaff) == "3.4 ax + 5.1"
        @test string(2 - aaff) == "-3.4 ax - 1.1"
        @test string(2 * aaff) == "6.8 ax + 6.2"
        @test_throws ErrorException 2 / aaff
    end
end


@testset "Variable--... tests" begin
    # Variable--Uncertain
    @test string(x + a) == "x + a"
    @test string(x - a) == "x + -a"
    @test string(x * a) == "a x"
    @test_throws ErrorException string(x / a)
    # Variable--UncExpr
    @test string(x + uaff) == "x + 2.3 a + 5.5"
    @test string(x - uaff) == "x + -2.3 a - 5.5"
    @test string(x * uaff) == "(2.3 a + 5.5) x"
    @test_throws ErrorException string(x / uaff)
    # Variable--UncVarExpr
    @test string(x + faff) == "(5 a + 1) x + x + 2 b + 3"
    @test string(x - faff) == "(-5 a - 1) x + x + -2 b - 3"
    @test_throws ErrorException x * faff
    @test_throws ErrorException x / faff
    @testset "Adaptive" begin
        @test string(x + ax) == "x + ax"
        @test string(x - ax) == "x - ax"
        @test_throws ErrorException x * ax
        @test_throws ErrorException x / ax
    end
    @testset "AdaptExpr" begin
        @test string(x + aaff) == "x + 3.4 ax + 3.1"
        @test string(x - aaff) == "x - 3.4 ax - 3.1"
        @test_throws ErrorException x * aaff
        @test_throws ErrorException x / aaff
    end
end


@testset "AffExpr--... tests" begin  # aff  = 7.1 * x + 2.5
    # AffExpr--Uncertain
    @test string(aff + a) == "7.1 x + a + 2.5"
    @test string(aff - a) == "7.1 x + -a + 2.5"
    @test string(aff * a) == "(7.1 a) x + 2.5 a"
    @test_throws ErrorException aff / a
    # AffExpr--UncExpr
    @test string(aff + uaff) == "7.1 x + 2.3 a + 8"
    @test string(aff - uaff) == "7.1 x + -2.3 a - 3"
    @test string(aff * uaff) == "(16.33 a + 39.05) x + 5.75 a + 13.75"
    @test_throws ErrorException aff / uaff
    # AffExpr--UncVarExpr
    @test string(aff + faff) == "7.1 x + (5 a + 1) x + 2 b + 5.5"
    @test string(aff - faff) == "7.1 x + (-5 a - 1) x + -2 b - 0.5"
    @test_throws ErrorException aff * faff
    @test_throws ErrorException aff / faff
    @testset "Adaptive" begin
        @test string(aff + ax) == "7.1 x + ax + 2.5"
        @test string(aff - ax) == "7.1 x - ax + 2.5"
        @test_throws ErrorException aff * ax
        @test_throws ErrorException aff / ax
    end
    @testset "AdaptExpr" begin  # aaff = 3.4 * ax + 3.1
        @test string(aff + aaff) == "7.1 x + 3.4 ax + 5.6"
        @test string(aff - aaff) == "7.1 x - 3.4 ax - 0.6000000000000001"
        @test_throws ErrorException aff * aaff
        @test_throws ErrorException aff / aaff
    end
end


@testset "Uncertain--... tests" begin
    @test !isequal(a,b)
    # Uncertain--Number
    @test string(a + 4.13) == "a + 4.13"
    @test string(a - 3.16) == "a - 3.16"
    @test string(a * 5.23) == "5.23 a"
    @test string(a / 2.0) == "0.5 a"
    # Uncertain--Variable
    @test string(a + x) == "x + a"
    @test string(a - x) == "-x + a"
    @test string(a * x) == "a x"
    @test_throws ErrorException string(a / x)
    # Uncertain--AffExpr
    @test string(a + aff) == "7.1 x + a + 2.5"
    @test string(a - aff) == "-7.1 x + a - 2.5"
    @test string(a * aff) == "(7.1 a) x + 2.5 a"
    @test_throws ErrorException a / aff
    # Uncertain--Uncertain
    @test string(a + b) == "a + b"
    @test string(a - b) == "a - b"
    @test_throws ErrorException a * b
    @test_throws ErrorException a / b
    # Uncertain--UncExpr (uaff = 2.3 * a + 5.5)
    @test string(b + uaff) == "b + 2.3 a + 5.5"
    @test string(b - uaff) == "b - 2.3 a - 5.5"
    @test_throws ErrorException b * uaff
    @test_throws ErrorException b / uaff
    # Uncertain--UncVarExpr (faff = (5a + 1)x + 2b + 3)
    @test string(a + faff) == "(5 a + 1) x + a + 2 b + 3"
    @test string(a - faff) == "(-5 a - 1) x + a - 2 b - 3"
    @test_throws ErrorException a * faff
    @test_throws ErrorException a / faff
    @testset "Adaptive" begin
        @test string(b + ax) ==  "ax + b"
        @test string(b - ax) == "-ax + b"
        @test string(b * ax) == "b ax"
        @test_throws ErrorException b / ax
    end
    @testset "AdaptExpr" begin  # aaff = 3.4 * ax + 3.1
        @test string(b + aaff) ==  "3.4 ax + b + 3.1"
        @test string(b - aaff) == "-3.4 ax + b - 3.1"
        @test string(b * aaff) == "(3.4 b) ax + 3.1 b"
        @test_throws ErrorException b / aaff
    end
end  # "Uncertain--... tests"


@testset "UncExpr--... tests" begin
    # Constructors
    @test string(JuMPeR.UncExpr(1)) == "1"
    @test string(JuMPeR.UncExpr(a)) == "a"
    @test typeof(zero(uaff)) == JuMPeR.UncExpr
    @test string(JuMPeR.UncExpr(0)) == "0"
    # UncExpr--Number
    @test string(uaff + 4.0) == "2.3 a + 9.5"
    @test string(uaff - 3.0) == "2.3 a + 2.5"
    @test string(uaff * 2.0) == "4.6 a + 11"
    @test string(uaff / 2.0) == "1.15 a + 2.75"
    # UncExpr--Variable
    @test string(uaff + x) == "x + 2.3 a + 5.5"
    @test string(uaff - x) == "-x + 2.3 a + 5.5"
    @test string(uaff * x) == "(2.3 a + 5.5) x"
    @test_throws ErrorException uaff / x
    # UncExpr--AffExpr (aff = 7.1 x + 2.5)
    @test string(uaff + aff) == "7.1 x + 2.3 a + 8"
    @test string(uaff - aff) == "-7.1 x + 2.3 a + 3"
    @test string(uaff * aff) == "(16.33 a + 39.05) x + 5.75 a + 13.75"
    @test_throws ErrorException uaff / aff
    # UncExpr--Uncertain
    @test string(uaff + b) == "b + 2.3 a + 5.5"
    @test string(uaff - b) == "-b + 2.3 a + 5.5"
    @test_throws ErrorException uaff * b
    @test_throws ErrorException uaff / b
    # UncExpr--UncExpr (uaff2 = 3.4 b + 1.1)
    @test string(uaff + uaff2) == "2.3 a + 3.4 b + 6.6"
    @test string(uaff - uaff2) == "2.3 a - 3.4 b + 4.4"
    @test_throws ErrorException uaff * uaff2
    @test_throws ErrorException uaff / uaff2
    # UncExpr--UncVarExpr (faff = (5a + 1)x + 2b + 3)
    @test string(uaff + faff) == "(5 a + 1) x + 2.3 a + 2 b + 8.5"
    @test string(uaff - faff) == "(-5 a - 1) x + 2.3 a - 2 b + 2.5"
    @test_throws ErrorException uaff * faff
    @test_throws ErrorException uaff / faff
    @testset "Adaptive" begin  # uex  = 2.3 * b + 5.5
        @test string(uex + ax) ==  "ax + 2.3 b + 5.5"
        @test string(uex - ax) == "-ax + 2.3 b + 5.5"
        @test string(uex * ax) == "(2.3 b + 5.5) ax"
        @test_throws ErrorException uex / ax
    end
    @testset "AdaptExpr" begin  # aaff = 3.4 * ax + 3.1
        @test string(uex + aaff) ==  "3.4 ax + 2.3 b + 8.6"
        @test string(uex - aaff) == "-3.4 ax + 2.3 b + 2.4"
        @test string(uex * aaff) == "(7.819999999999999 b + 18.7) ax + 7.13 b + 17.05"
        @test_throws ErrorException uex / aaff
    end
end  # "UncExpr--... tests"


@testset "UncVarExpr--... tests" begin
    # Constructors
    @test typeof(zero(faff)) == JuMPeR.UncVarExpr
    @test string(JuMPeR.UncVarExpr()) == "0"
    # Push/append
    pusher = a * x
    @test string(pusher) == "a x"
    push!(pusher, 2.0, y)
    @test string(pusher) == "a x + 2 y"
    push!(pusher, uaff, z)
    @test string(pusher) == "a x + 2 y + (2.3 a + 5.5) z"
    # faff = (5a + 1)x + 2b + 3)
    # UncVarExpr--Number
    @test string(faff + 4.0) == "(5 a + 1) x + 2 b + 7"
    @test string(faff - 2.0) == "(5 a + 1) x + 2 b + 1"
    @test string(faff * 2.0) == "(10 a + 2) x + 4 b + 6"
    @test string(faff / 2.0) == "(2.5 a + 0.5) x + b + 1.5"
    # UncVarExpr--Variable
    @test string(faff + y) == "(5 a + 1) x + y + 2 b + 3"
    @test string(faff - y) == "(5 a + 1) x - y + 2 b + 3"
    @test_throws ErrorException faff * y
    @test_throws ErrorException faff / y
    # UncVarExpr--AffExpr (aff2 = 1.2y + 1.2)
    @test string(faff + aff2) == "1.2 y + (5 a + 1) x + 2 b + 4.2"
    @test string(faff - aff2) == "(5 a + 1) x - 1.2 y + 2 b + 1.8"
    @test_throws ErrorException faff * aff2
    @test_throws ErrorException faff / aff2
    # UncVarExpr--Uncertain
    @test string(faff + a) == "(5 a + 1) x + a + 2 b + 3"
    @test string(faff - a) == "(5 a + 1) x + -a + 2 b + 3"
    @test_throws ErrorException faff * a
    @test_throws ErrorException faff / a
    # UncVarExpr--UncExpr (uaff = 2.3 * a + 5.5)
    @test string(faff + uaff) == "(5 a + 1) x + 2 b + 2.3 a + 8.5"
    @test string(faff - uaff) == "(5 a + 1) x + 2 b - 2.3 a - 2.5"
    @test_throws ErrorException faff * uaff
    @test_throws ErrorException faff / uaff
    # UncVarExpr--UncVarExpr
    @test string(faff + faff) == "(5 a + 1) x + (5 a + 1) x + 4 b + 6"
    @test string(faff - faff) == "(5 a + 1) x + (-5 a - 1) x"
    @test_throws ErrorException faff * faff
    @test_throws ErrorException faff / faff
    @testset "Adaptive" begin  # uaex = (5b+1)x + (2c + 3)
        @test string(uaex + ax) == "(5 b + 1) x + ax + 2 c + 3"
        @test string(uaex - ax) == "(5 b + 1) x - ax + 2 c + 3"
        @test_throws ErrorException uaex * ax
        @test_throws ErrorException uaex / ax
    end
    @testset "AdaptExpr" begin  # aaff = 3.4 * ax + 3.1
        @test string(uaex + aaff) == "(5 b + 1) x + 3.4 ax + 2 c + 6.1"
        @test string(uaex - aaff) == "(5 b + 1) x - 3.4 ax + 2 c - 0.10000000000000009"
        @test_throws ErrorException uaex * aaff
        @test_throws ErrorException uaex / aaff
    end
end  # "UncVarExpr--... tests"


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

    @testset "UncVarExpr" begin
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
    @testset "AdaptExpr" begin
        @test string(ay + aaff) == "ay + 3.4 ax + 3.1"
        @test string(ay - aaff) == "ay - 3.4 ax - 3.1"
        @test_throws ErrorException ay * aaff
        @test_throws ErrorException ay / aaff
    end
end  # "Adaptive--"


@testset "AdaptExpr--" begin
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
    @testset "UncVarExpr" begin
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
end  # "AdaptExpr--"

end  # "Basics"


@testset "Higher level" begin

m = RobustModel()
@variable(m, 0 <= x[1:3] <= 1)
@uncertain(m, 2 <= a <= 3)
@uncertain(m, 5 <= b <= 6)
@uncertain(m, u[1:3])
@uncertain(m, v[4:6])
@uncertain(m, U[1:3,1:3])
@uncertain(m, w[[:foo,:bar]])

@testset "sum" begin
    @test string(sum(u)) == "u[3] + u[1] + u[2]"
    @test string(sum(U)) == "U[3,3] + U[2,3] + U[1,3] + U[3,2] + U[2,2] + U[1,2] + U[3,1] + U[1,1] + U[2,1]"
    @test string(sum(w)) in ["w[foo] + w[bar]", "w[bar] + w[foo]"]
    @test string(sum([2.0*a, 4.0*b])) == "2 a + 4 b"
    @test string(sum([x[1] + 2.0*a, x[2] + 4.0*b])) == "x[1] + x[2] + 2 a + 4 b"
end

@testset "dot" begin
    c = [3.5, 4.0, 2.0]
    A = [3.0 4.0 5.0;
         1.5 2.5 3.3;
         5.5 6.2 1.2]
    # DOT
    # Vector{Float64} :: JuMP.JuMPArray{Uncertain}
    @test string(dot(c,u)) == "3.5 u[1] + 4 u[2] + 2 u[3]"
    @test string(dot(u,c)) == "3.5 u[1] + 4 u[2] + 2 u[3]"
    # Matrix{Float64} (2D) :: JuMP.JuMPArray{Uncertain} (2D)
    @test string(dot(A,U)) == "3 U[1,1] + 1.5 U[2,1] + 5.5 U[3,1] + 4 U[1,2] + 2.5 U[2,2] + 6.2 U[3,2] + 5 U[1,3] + 3.3 U[2,3] + 1.2 U[3,3]"

    # JuMP.JuMPArray{Variable} :: JuMP.JuMPArray{Uncertain}
    @test string(dot(x,u)) == "u[1] x[1] + u[2] x[2] + u[3] x[3]"

    # Matrix{Float64} (2D) :: JuMP.JuMPDict{Uncertain} (1D)
    @test_throws DimensionMismatch dot(A, u)
end

end  # "Higher level"


@testset "Matrix operations" begin

    m = RobustModel()
    @variable(m, matvar[1:3,1:3])
    @uncertain(m, vecunc[1:3])
    @uncertain(m, matunc[1:3,1:3])
    s = JuMPeR.get_robust(m).default_uncset

    b = [1,2,3]
    A = Matrix(I,3,3)
    @constraint(m, A*vecunc .== b)
    c = s.linear_constraints[1:3]
    for i in 1:3
        @test string(c[i]) == "vecunc[$i] $eq $i"
    end

    @constraint(m, matunc*ones(3) .== b)
    c = s.linear_constraints[4:6]
    for i in 1:3
        @test string(c[i]) == "matunc[$i,1] + matunc[$i,2] + matunc[$i,3] $eq $i"
    end

    @constraint(m, matvar*vecunc .== b)
    c = JuMPeR.get_robust(m).unc_constraints[1:3]
    for i in 1:3
        @test string(c[i]) == "vecunc[1] matvar[$i,1] + vecunc[2] matvar[$i,2] + vecunc[3] matvar[$i,3] $eq $i"
    end

    @constraint(m, matvar*matunc .== ones(3,3))
    c = reshape(JuMPeR.get_robust(m).unc_constraints[4:12], (3,3))
    for i in 1:3, j in 1:3
        @test string(c[i,j]) == "matunc[1,$j] matvar[$i,1] + matunc[2,$j] matvar[$i,2] + matunc[3,$j] matvar[$i,3] $eq 1"
    end

    @constraint(m, matvar.*matunc .== ones(3,3))
    c = reshape(JuMPeR.get_robust(m).unc_constraints[13:21], (3,3))
    for i in 1:3, j in 1:3
        @test string(c[i,j]) == "matunc[$i,$j] matvar[$i,$j] $eq 1"
    end

    @constraint(m, 2 .* matvar.*matunc + matvar.*matunc .== ones(3,3))
    c = reshape(JuMPeR.get_robust(m).unc_constraints[22:30], (3,3))
    for i in 1:3, j in 1:3
        @test string(c[i,j]) == "(2 matunc[$i,$j]) matvar[$i,$j] + matunc[$i,$j] matvar[$i,$j] $eq 1"
    end
end  # "Matrix operations"


@testset "Unc. set norms" begin

rm = RobustModel()
us = JuMPeR.get_robust(rm).default_uncset
nc = us.norm_constraints
@uncertain(rm, u[1:3])
@constraint(us, norm((u[i] for i=1:3),1) <= 1)
@test string(nc[end]) == "‖u[1],u[2],u[3]‖₁ $leq 1"
@constraint(us, norm(u[i] for i=1:3) <= 2)
@test string(nc[end]) == "‖u[1],u[2],u[3]‖₂ $leq 2"
@constraint(us, norm((u[i] for i=1:3),Inf) <= 9)
@test string(nc[end]) == "‖u[1],u[2],u[3]‖∞ $leq 9"

@constraint(us, 2*norm((u[i] for i=1:3),1) <= 1)
@test string(nc[end]) == "2‖u[1],u[2],u[3]‖₁ $leq 1"
@constraint(us, -1*norm((u[i] for i=1:3),1) >= -1)
@test string(nc[end]) == "‖u[1],u[2],u[3]‖₁ $leq 1"

@constraint(us, 1 + norm((u[i] for i=1:3),1) <= 1)
@test string(nc[end]) == "‖u[1],u[2],u[3]‖₁ $leq 0"

@variable(rm, x)

@test_throws ErrorException @constraint(us, norm((u[i] for i=1:3),1) + u[1] <= 1)
@test_throws ErrorException @constraint(us, norm((u[i] for i=1:3),1) - u[1] <= 1)
@test_throws ErrorException @constraint(us, norm((u[i] for i=1:3),1) * u[1] <= 1)
@test_throws ErrorException @constraint(us, norm((u[i] for i=1:3),1) / u[1] <= 1)

@test_throws ErrorException @constraint(us, norm((u[i] for i=1:3),1) + (2*u[1]) <= 1)
@test_throws ErrorException @constraint(us, norm((u[i] for i=1:3),1) - (2*u[1]) <= 1)
@test_throws ErrorException @constraint(us, norm((u[i] for i=1:3),1) * (2*u[1]) <= 1)
@test_throws MethodError    @constraint(us, norm((u[i] for i=1:3),1) / (2*u[1]) <= 1)
# MethodError: `/` has no method matching /(::Int64, ::JuMP.GenericAffExpr{Float64,JuMPeR.Uncertain})

@test_throws ErrorException @constraint(us, norm((u[i] for i=1:3),1) + (2*u[1]*x+u[2]) <= 1)
@test_throws ErrorException @constraint(us, norm((u[i] for i=1:3),1) - (2*u[1]*x+u[2]) <= 1)
@test_throws ErrorException @constraint(us, norm((u[i] for i=1:3),1) * (2*u[1]*x+u[2]) <= 1)
@test_throws MethodError    @constraint(us, norm((u[i] for i=1:3),1) / (2*u[1]*x+u[2]) <= 1)

@test_throws ErrorException @constraint(us, x + norm((u[i] for i=1:3),1) <= 1)
@test_throws ErrorException @constraint(us, x - norm((u[i] for i=1:3),1) <= 1)
@test_throws ErrorException @constraint(us, (2x) + norm((u[i] for i=1:3),1) <= 1)
@test_throws ErrorException @constraint(us, (2x) - norm((u[i] for i=1:3),1) <= 1)

@test_throws ErrorException @constraint(us, (u[1]) + norm((u[i] for i=1:3),1) <= 1)
@test_throws ErrorException @constraint(us, (u[1]) - norm((u[i] for i=1:3),1) <= 1)
@test_throws ErrorException @constraint(us, (u[1]) * norm((u[i] for i=1:3),1) <= 1)
# @constraint(us, (u[1]) / norm((u[i] for i=1:3),1) <= 1)
# UndefVarError: i not defined

@test_throws ErrorException @constraint(us, (2*u[1]) + norm((u[i] for i=1:3),1) <= 1)
@test_throws ErrorException @constraint(us, (2*u[1]) - norm((u[i] for i=1:3),1) <= 1)
@test_throws ErrorException @constraint(us, (2*u[1]) * norm((u[i] for i=1:3),1) <= 1)
# @constraint(us, (2*u[1]) / norm((u[i] for i=1:3),1) <= 1)
# UndefVarError: i not defined

@test_throws ErrorException @constraint(us, (2*u[1]*x+u[2]) + norm((u[i] for i=1:3),1) <= 1)
@test_throws ErrorException @constraint(us, (2*u[1]*x+u[2]) - norm((u[i] for i=1:3),1) <= 1)
@test_throws ErrorException @constraint(us, (2*u[1]*x+u[2]) * norm((u[i] for i=1:3),1) <= 1)
# @constraint(us, (2*u[1]*x+u[2]) / norm((u[i] for i=1:3),1) <= 1)
# UndefVarError: i not defined

end  # "Unc. set norms"


end  # "Operators"
