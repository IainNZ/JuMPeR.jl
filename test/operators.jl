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
using BaseTestNext
import JuMP: REPLMode, IJuliaMode

# To ensure the tests work on both Windows and Linux/OSX, we need
# to use the correct comparison operators in the strings
const leq = JuMP.repl[:leq]
const geq = JuMP.repl[:geq]
const  eq = JuMP.repl[:eq]

@testset "Operators" begin
print_with_color(:yellow, "Operators...\n")

@testset "Basics" begin

# Setup
m = RobustModel()
@defVar(m, w)
@defVar(m, x)
@defVar(m, y)
@defVar(m, z)

@defUnc(m, a)
@defUnc(m, b)

@test getNumUncs(m) == 2
@test_throws ErrorException getNumUncs(Model())

aff   = 7.1 * x + 2.5
aff2  = 1.2 * y + 1.2
uaff  = 2.3 * a + 5.5
uaff2 = 3.4 * b + 1.1
faff  = (5a+1)x + (2b + 3)

@test sprint(print, a)  == "a"
@test sprint( show, a)  == "a"
@test affToStr(aff)     == "7.1 x + 2.5"
@test affToStr(aff2)    == "1.2 y + 1.2"
@test affToStr(uaff)    == "2.3 a + 5.5"
@test affToStr(uaff2)   == "3.4 b + 1.1"
@test affToStr(faff)    == "(5 a + 1) x + 2 b + 3"
@test affToStr(x+UncExpr(0.0)) == sprint(print, x)

@testset "Number--... tests" begin
    # Number--Uncertain
    @test affToStr(4.13 + a) == "a + 4.13"
    @test affToStr(3.16 - a) == "-a + 3.16"
    @test affToStr(5.23 * a) == "5.23 a"
    @test_throws ErrorException 2.94 / a
    # Number--UncExpr
    @test affToStr(2.3 + uaff) == "2.3 a + 7.8"
    @test affToStr(1.5 - uaff) == "-2.3 a - 4"
    @test affToStr(2.0 * uaff) == "4.6 a + 11"
    @test_throws MethodError 2.94 / uaff
    # Number--UncAffExpr
    @test affToStr(2.3 + faff) == "(5 a + 1) x + 2 b + 5.3"
    @test affToStr(1.0 - faff) == "(-5 a - 1) x + -2 b - 2"
    @test affToStr(2.0 * faff) == "(10 a + 2) x + 4 b + 6"
    @test_throws MethodError 2.94 / faff
end


@testset "Variable--... tests" begin
    # Variable--Uncertain
    @test affToStr(x + a) == "x + a"
    @test affToStr(x - a) == "x + -a"
    @test affToStr(x * a) == "a x"
    @test_throws ErrorException affToStr(x / a)
    # Variable--UncExpr
    @test affToStr(x + uaff) == "x + 2.3 a + 5.5"
    @test affToStr(x - uaff) == "x + -2.3 a - 5.5"
    @test affToStr(x * uaff) == "(2.3 a + 5.5) x"
    @test_throws ErrorException affToStr(x / uaff)
    # Variable--UncAffExpr
    @test affToStr(x + faff) == "(5 a + 1) x + x + 2 b + 3"
    @test affToStr(x - faff) == "(-5 a - 1) x + x + -2 b - 3"
    @test_throws ErrorException x * faff
    @test_throws ErrorException x / faff
end


@testset "AffExpr--... tests" begin
    # AffExpr--Uncertain
    @test affToStr(aff + a) == "7.1 x + a + 2.5"
    @test affToStr(aff - a) == "7.1 x + -a + 2.5"
    @test affToStr(aff * a) == "(7.1 a) x + 2.5 a"
    @test_throws ErrorException aff / a
    # AffExpr--UncExpr
    @test affToStr(aff + uaff) == "7.1 x + 2.3 a + 8"
    @test affToStr(aff - uaff) == "7.1 x + -2.3 a - 3"
    @test affToStr(aff * uaff) == "(16.33 a + 39.05) x + 5.75 a + 13.75"
    @test_throws ErrorException aff / uaff
    # AffExpr--UncAffExpr
    @test affToStr(aff + faff) == "7.1 x + (5 a + 1) x + 2 b + 5.5"
    @test affToStr(aff - faff) == "7.1 x + (-5 a - 1) x + -2 b - 0.5"
    @test_throws ErrorException aff * faff
    @test_throws ErrorException aff / faff
end


@testset "Uncertain--... tests" begin
    @test !isequal(a,b)
    @test getName(a) == "a"
    # Uncertain--Number
    @test affToStr(a + 4.13) == "a + 4.13"
    @test affToStr(a - 3.16) == "a - 3.16"
    @test affToStr(a * 5.23) == "5.23 a"
    @test affToStr(a / 2.0) == "0.5 a"
    # Uncertain--Variable
    @test affToStr(a + x) == "x + a"
    @test affToStr(a - x) == "-x + a"
    @test affToStr(a * x) == "a x"
    @test_throws ErrorException affToStr(a / x)
    # Uncertain--AffExpr
    @test affToStr(a + aff) == "7.1 x + a + 2.5"
    @test affToStr(a - aff) == "-7.1 x + a - 2.5"
    @test affToStr(a * aff) == "(7.1 a) x + 2.5 a"
    @test_throws ErrorException a / aff
    # Uncertain--Uncertain
    @test affToStr(a + b) == "a + b"
    @test affToStr(a - b) == "a - b"
    @test_throws ErrorException a * b
    @test_throws ErrorException a / b
    # Uncertain--UncExpr (uaff = 2.3 * a + 5.5)
    @test affToStr(b + uaff) == "b + 2.3 a + 5.5"
    @test affToStr(b - uaff) == "b - 2.3 a - 5.5"
    @test_throws ErrorException b * uaff
    @test_throws ErrorException b / uaff
    # Uncertain--UncAffExpr (faff = (5a + 1)x + 2b + 3)
    @test affToStr(a + faff) == "(5 a + 1) x + a + 2 b + 3"
    @test affToStr(a - faff) == "(-5 a - 1) x + a - 2 b - 3"
    @test_throws ErrorException a * faff
    @test_throws ErrorException a / faff
end  # "Uncertain--... tests"


@testset "UncExpr--... tests" begin
    # Constructors
    @test affToStr(UncExpr()) == "0"
    @test affToStr(UncExpr(1)) == "1"
    @test affToStr(UncExpr(a)) == "a"
    @test typeof(zero(uaff)) == UncExpr
    @test affToStr(zero(UncExpr)) == "0"
    # UncExpr--Number
    @test affToStr(uaff + 4.0) == "2.3 a + 9.5"
    @test affToStr(uaff - 3.0) == "2.3 a + 2.5"
    @test affToStr(uaff * 2.0) == "4.6 a + 11"
    @test affToStr(uaff / 2.0) == "1.15 a + 2.75"
    # UncExpr--Variable
    @test affToStr(uaff + x) == "x + 2.3 a + 5.5"
    @test affToStr(uaff - x) == "-x + 2.3 a + 5.5"
    @test affToStr(uaff * x) == "(2.3 a + 5.5) x"
    @test_throws ErrorException uaff / x
    # UncExpr--AffExpr (aff = 7.1 x + 2.5)
    @test affToStr(uaff + aff) == "7.1 x + 2.3 a + 8"
    @test affToStr(uaff - aff) == "-7.1 x + 2.3 a + 3"
    @test affToStr(uaff * aff) == "(16.33 a + 39.05) x + 5.75 a + 13.75"
    @test_throws ErrorException uaff / aff
    # UncExpr--Uncertain
    @test affToStr(uaff + b) == "b + 2.3 a + 5.5"
    @test affToStr(uaff - b) == "-b + 2.3 a + 5.5"
    @test_throws ErrorException uaff * b
    @test_throws ErrorException uaff / b
    # UncExpr--UncExpr (uaff2 = 3.4 b + 1.1)
    @test affToStr(uaff + uaff2) == "2.3 a + 3.4 b + 6.6"
    @test affToStr(uaff - uaff2) == "2.3 a - 3.4 b + 4.4"
    @test_throws ErrorException uaff * uaff2
    @test_throws ErrorException uaff / uaff2
    # UncExpr--UncAffExpr (faff = (5a + 1)x + 2b + 3)
    @test affToStr(uaff + faff) == "(5 a + 1) x + 2.3 a + 2 b + 8.5"
    @test affToStr(uaff - faff) == "(-5 a - 1) x + 2.3 a - 2 b + 2.5"
    @test_throws ErrorException uaff * faff
    @test_throws ErrorException uaff / faff
end  # "UncExpr--... tests"


@testset "UncAffExpr--... tests" begin
    # Constructors
    @test affToStr(UncAffExpr()) == "0"
    @test typeof(zero(faff)) == UncAffExpr
    @test affToStr(zero(UncAffExpr)) == "0"
    # Push/append
    pusher = a * x
    @test affToStr(pusher) == "a x"
    push!(pusher, 2.0, y)
    @test affToStr(pusher) == "a x + 2 y"
    push!(pusher, uaff, z)
    @test affToStr(pusher) == "a x + 2 y + (2.3 a + 5.5) z"
    # faff = (5a + 1)x + 2b + 3)
    # UncAffExpr--Number
    @test affToStr(faff + 4.0) == "(5 a + 1) x + 2 b + 7"
    @test affToStr(faff - 2.0) == "(5 a + 1) x + 2 b + 1"
    @test affToStr(faff * 2.0) == "(10 a + 2) x + 4 b + 6"
    @test affToStr(faff / 2.0) == "(2.5 a + 0.5) x + b + 1.5"
    # UncAffExpr--Variable
    @test affToStr(faff + y) == "(5 a + 1) x + y + 2 b + 3"
    @test affToStr(faff - y) == "(5 a + 1) x - y + 2 b + 3"
    @test_throws ErrorException faff * y
    @test_throws ErrorException faff / y
    # UncAffExpr--AffExpr (aff2 = 1.2y + 1.2)
    @test affToStr(faff + aff2) == "1.2 y + (5 a + 1) x + 2 b + 4.2"
    @test affToStr(faff - aff2) == "(5 a + 1) x - 1.2 y + 2 b + 1.8"
    @test_throws ErrorException faff * aff2
    @test_throws ErrorException faff / aff2
    # UncAffExpr--Uncertain
    @test affToStr(faff + a) == "(5 a + 1) x + a + 2 b + 3"
    @test affToStr(faff - a) == "(5 a + 1) x + -a + 2 b + 3"
    @test_throws ErrorException faff * a
    @test_throws ErrorException faff / a
    # UncAffExpr--UncExpr (uaff = 2.3 * a + 5.5)
    @test affToStr(faff + uaff) == "(5 a + 1) x + 2 b + 2.3 a + 8.5"
    @test affToStr(faff - uaff) == "(5 a + 1) x + 2 b - 2.3 a - 2.5"
    @test_throws ErrorException faff * uaff
    @test_throws ErrorException faff / uaff
    # UncAffExpr--UncAffExpr
    @test affToStr(faff + faff) == "(5 a + 1) x + (5 a + 1) x + 4 b + 6"
    @test affToStr(faff - faff) == "(5 a + 1) x + (-5 a - 1) x"
    @test_throws ErrorException faff * faff
    @test_throws ErrorException faff / faff
end  # "UncAffExpr--... tests"

end  # "Basics"


@testset "Higher level" begin

m = RobustModel()
@defVar(m, 0 <= x[1:3] <= 1)
@defUnc(m, 2 <= a <= 3)
@defUnc(m, 5 <= b <= 6)
@defUnc(m, u[1:3])
@defUnc(m, v[4:6])
@defUnc(m, U[1:3,1:3])
@defUnc(m, w[[:foo,:bar]])

@testset "sum" begin
    @test affToStr(sum(u)) == "u[3] + u[1] + u[2]"
    @test affToStr(sum(U)) == "U[3,3] + U[2,3] + U[1,3] + U[3,2] + U[2,2] + U[1,2] + U[3,1] + U[1,1] + U[2,1]"
    @test affToStr(sum(w)) in ["w[foo] + w[bar]", "w[bar] + w[foo]"]
    @test affToStr(sum([2.0*a, 4.0*b])) == "2 a + 4 b"
    @test affToStr(sum([x[1] + 2.0*a, x[2] + 4.0*b])) == "x[1] + x[2] + 2 a + 4 b"
end

@testset "dot" begin
    c = [3.5, 4.0, 2.0]
    A = [3.0 4.0 5.0;
         1.5 2.5 3.3;
         5.5 6.2 1.2]
    # DOT
    # Vector{Float64} :: JuMPArray{Uncertain}
    @test affToStr(dot(c,u)) == "3.5 u[1] + 4 u[2] + 2 u[3]"
    @test affToStr(dot(u,c)) == "3.5 u[1] + 4 u[2] + 2 u[3]"
    # Matrix{Float64} (2D) :: JuMPArray{Uncertain} (2D)
    @test affToStr(vecdot(A,U)) == "3 U[1,1] + 1.5 U[2,1] + 5.5 U[3,1] + 4 U[1,2] + 2.5 U[2,2] + 6.2 U[3,2] + 5 U[1,3] + 3.3 U[2,3] + 1.2 U[3,3]"

    # JuMPArray{Variable} :: JuMPArray{Uncertain}
    @test affToStr(dot(x,u)) == "u[1] x[1] + u[2] x[2] + u[3] x[3]"

    # Matrix{Float64} (2D) :: JuMPDict{Uncertain} (1D)
    @test_throws DimensionMismatch vecdot(A, u)
end

end  # "Higher level"


@testset "Matrix operations" begin

    m = RobustModel()
    @defVar(m, matvar[1:3,1:3])
    @defUnc(m, vecunc[1:3])
    @defUnc(m, matunc[1:3,1:3])

    b = [1,2,3]
    A = eye(3,3)
    c = @addConstraint(m, A*vecunc .== b)
    for i in 1:3
        @test conToStr(c[i]) == "vecunc[$i] $eq $i"
    end

    c = @addConstraint(m, matunc*ones(3) .== b)
    for i in 1:3
        @test conToStr(c[i]) == "matunc[$i,1] + matunc[$i,2] + matunc[$i,3] $eq $i"
    end

    c = @addConstraint(m, matvar*vecunc .== b)
    for i in 1:3
        @test conToStr(c[i]) == "vecunc[1] matvar[$i,1] + vecunc[2] matvar[$i,2] + vecunc[3] matvar[$i,3] $eq $i"
    end

    c = @addConstraint(m, matvar*matunc .== ones(3,3))
    for i in 1:3, j in 1:3
        @test conToStr(c[i,j]) == "matunc[1,$j] matvar[$i,1] + matunc[2,$j] matvar[$i,2] + matunc[3,$j] matvar[$i,3] $eq 1"
    end

    c = @addConstraint(m, matvar.*matunc .== ones(3,3))
    for i in 1:3, j in 1:3
        @test conToStr(c[i,j]) == "matunc[$i,$j] matvar[$i,$j] $eq 1"
    end

    c = @addConstraint(m, 2.*matvar.*matunc + matvar.*matunc .== ones(3,3))
    for i in 1:3, j in 1:3
        @test conToStr(c[i,j]) == "(2 matunc[$i,$j]) matvar[$i,$j] + matunc[$i,$j] matvar[$i,$j] $eq 1"
    end
end  # "Matrix operations"


@testset "Unc. set norms" begin

rm = RobustModel()
nc = JuMPeR.getRobust(rm).normconstraints
@defUnc(rm, u[1:3])
@addConstraint(rm, norm1{u[i],i=1:3} <= 1)
@test conToStr(nc[end]) == "‖u[1],u[2],u[3]‖₁ $leq 1"
@addConstraint(rm, norm2{u[i],i=1:3} <= 2)
@test conToStr(nc[end]) == "‖u[1],u[2],u[3]‖₂ $leq 2"
@addConstraint(rm, norm∞{u[i],i=1:3} <= 9)
@test conToStr(nc[end]) == "‖u[1],u[2],u[3]‖∞ $leq 9"

@addConstraint(rm, 2*norm1{u[i],i=1:3} <= 1)
@test conToStr(nc[end]) == "2‖u[1],u[2],u[3]‖₁ $leq 1"
@addConstraint(rm, -1*norm1{u[i],i=1:3} >= -1)
@test conToStr(nc[end]) == "‖u[1],u[2],u[3]‖₁ $leq 1"

@addConstraint(rm, 1 + norm1{u[i],i=1:3} <= 1)
@test conToStr(nc[end]) == "‖u[1],u[2],u[3]‖₁ $leq 0"

@defVar(rm, x)

@test_throws ErrorException @addConstraint(rm, norm1{u[i],i=1:3} + u[1] <= 1)
@test_throws ErrorException @addConstraint(rm, norm1{u[i],i=1:3} - u[1] <= 1)
@test_throws ErrorException @addConstraint(rm, norm1{u[i],i=1:3} * u[1] <= 1)
@test_throws ErrorException @addConstraint(rm, norm1{u[i],i=1:3} / u[1] <= 1)

@test_throws ErrorException @addConstraint(rm, norm1{u[i],i=1:3} + (2*u[1]) <= 1)
@test_throws ErrorException @addConstraint(rm, norm1{u[i],i=1:3} - (2*u[1]) <= 1)
@test_throws ErrorException @addConstraint(rm, norm1{u[i],i=1:3} * (2*u[1]) <= 1)
@test_throws MethodError    @addConstraint(rm, norm1{u[i],i=1:3} / (2*u[1]) <= 1)
# MethodError: `/` has no method matching /(::Int64, ::JuMP.GenericAffExpr{Float64,JuMPeR.Uncertain})

@test_throws ErrorException @addConstraint(rm, norm1{u[i],i=1:3} + (2*u[1]*x+u[2]) <= 1)
@test_throws ErrorException @addConstraint(rm, norm1{u[i],i=1:3} - (2*u[1]*x+u[2]) <= 1)
@test_throws ErrorException @addConstraint(rm, norm1{u[i],i=1:3} * (2*u[1]*x+u[2]) <= 1)
@test_throws MethodError    @addConstraint(rm, norm1{u[i],i=1:3} / (2*u[1]*x+u[2]) <= 1)

@test_throws ErrorException @addConstraint(rm, x + norm1{u[i],i=1:3} <= 1)
@test_throws ErrorException @addConstraint(rm, x - norm1{u[i],i=1:3} <= 1)
@test_throws ErrorException @addConstraint(rm, (2x) + norm1{u[i],i=1:3} <= 1)
@test_throws ErrorException @addConstraint(rm, (2x) - norm1{u[i],i=1:3} <= 1)

@test_throws ErrorException @addConstraint(rm, (u[1]) + norm1{u[i],i=1:3} <= 1)
@test_throws ErrorException @addConstraint(rm, (u[1]) - norm1{u[i],i=1:3} <= 1)
@test_throws ErrorException @addConstraint(rm, (u[1]) * norm1{u[i],i=1:3} <= 1)
# @addConstraint(rm, (u[1]) / norm1{u[i],i=1:3} <= 1)
# UndefVarError: i not defined

@test_throws ErrorException @addConstraint(rm, (2*u[1]) + norm1{u[i],i=1:3} <= 1)
@test_throws ErrorException @addConstraint(rm, (2*u[1]) - norm1{u[i],i=1:3} <= 1)
@test_throws ErrorException @addConstraint(rm, (2*u[1]) * norm1{u[i],i=1:3} <= 1)
# @addConstraint(rm, (2*u[1]) / norm1{u[i],i=1:3} <= 1)
# UndefVarError: i not defined

@test_throws ErrorException @addConstraint(rm, (2*u[1]*x+u[2]) + norm1{u[i],i=1:3} <= 1)
@test_throws ErrorException @addConstraint(rm, (2*u[1]*x+u[2]) - norm1{u[i],i=1:3} <= 1)
@test_throws ErrorException @addConstraint(rm, (2*u[1]*x+u[2]) * norm1{u[i],i=1:3} <= 1)
# @addConstraint(rm, (2*u[1]*x+u[2]) / norm1{u[i],i=1:3} <= 1)
# UndefVarError: i not defined

end  # "Unc. set norms"


end  # "Operators"
