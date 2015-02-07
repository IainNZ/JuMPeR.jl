#############################################################################
# JuMPeR
# Julia for Mathematical Programming - extension for Robust Optimization
# See http://github.com/IainNZ/JuMPeR.jl
#############################################################################
# test/operators.jl
# Testing for all operator overloads
#############################################################################

using JuMP, JuMPeR
using FactCheck

# To ensure the tests work on both Windows and Linux/OSX, we need
# to use the correct comparison operators in the strings
const leq = JuMP.repl_leq
const geq = JuMP.repl_geq
const  eq = JuMP.repl_eq

facts("[operators] Robust operator tests") do

# Setup
m = RobustModel()
@defVar(m, w)
@defVar(m, x)
@defVar(m, y)
@defVar(m, z)

@defUnc(m, a)
@defUnc(m, b)

aff   = 7.1 * x + 2.5
aff2  = 1.2 * y + 1.2
uaff  = 2.3 * a + 5.5
uaff2 = 3.4 * b + 1.1
faff  = (5a+1)x + (2b + 3)

@fact affToStr(aff)   => "7.1 x + 2.5"
@fact affToStr(aff2)  => "1.2 y + 1.2"
@fact affToStr(uaff)  => "2.3 a + 5.5"
@fact affToStr(uaff2) => "3.4 b + 1.1"
@fact affToStr(faff)  => "(5 a + 1) x + 2 b + 3"

context("Core JuMPeR type methods") do
    # Model
    @fact getNumUncs(m) => 2
    @fact_throws getNumUncs(Model())

    # Uncertain
    @fact isequal(a,b) => false
    @fact getName(a) => "a"

    # UAffExpr
    @fact affToStr(UAffExpr()) => "0"
    @fact affToStr(UAffExpr(1)) => "1"
    @fact affToStr(UAffExpr(a)) => "a"
    @fact affToStr(UAffExpr(a,2)) => "2 a"
    @fact typeof(zero(uaff)) => UAffExpr
    @fact affToStr(zero(UAffExpr)) => "0"

    # FullAffExpr
    @fact affToStr(FullAffExpr()) => "0"
    @fact typeof(zero(faff)) => FullAffExpr
    @fact affToStr(zero(FullAffExpr)) => "0"
    pusher = a * x
    @fact affToStr(pusher) => "a x"
    push!(pusher, 2.0, y)
    @fact affToStr(pusher) => "a x + 2 y"
    push!(pusher, uaff, z)
    @fact affToStr(pusher) => "a x + 2 y + (2.3 a + 5.5) z"
end
    

context("Number--... tests") do
    # Number--Uncertain
    @fact affToStr(4.13 + a) => "a + 4.13"
    @fact affToStr(3.16 - a) => "-a + 3.16"
    @fact affToStr(5.23 * a) => "5.23 a"
    @fact_throws 2.94 / a
    # Number--UAffExpr
    @fact affToStr(2.3 + uaff) => "2.3 a + 7.8"
    @fact affToStr(1.5 - uaff) => "-2.3 a - 4"
    @fact affToStr(2.0 * uaff) => "4.6 a + 11"
    @fact_throws 2.94 / uaff
    # Number--FullAffExpr
    @fact affToStr(2.3 + faff) => "(5 a + 1) x + 2 b + 5.3"
    @fact affToStr(1.0 - faff) => "(-5 a - 1) x + -2 b - 2"
    @fact affToStr(2.0 * faff) => "(10 a + 2) x + 4 b + 6"
    @fact_throws 2.94 / faff
end


context("Variable--... tests") do
    # Variable--Uncertain
    @fact affToStr(x + a) => "x + a"
    @fact affToStr(x - a) => "x + -a"
    @fact affToStr(x * a) => "a x"
    @fact_throws affToStr(x / a)
    # Variable--UAffExpr
    @fact affToStr(x + uaff) => "x + 2.3 a + 5.5"
    @fact affToStr(x - uaff) => "x + -2.3 a - 5.5"
    @fact affToStr(x * uaff) => "(2.3 a + 5.5) x"
    @fact_throws affToStr(x / uaff)
    # Variable--FullAffExpr
    @fact affToStr(x + faff) => "(5 a + 1) x + x + 2 b + 3"
    @fact affToStr(x - faff) => "(-5 a - 1) x + x + -2 b - 3"
    @fact_throws x * faff
    @fact_throws x / faff
end


context("AffExpr--... tests") do
    # AffExpr--Uncertain
    @fact affToStr(aff + a) => "7.1 x + a + 2.5"
    @fact affToStr(aff - a) => "7.1 x + -a + 2.5"
    @fact affToStr(aff * a) => "(7.1 a) x + 2.5 a"
    @fact_throws aff / a
    # AffExpr--UAffExpr
    @fact affToStr(aff + uaff) => "7.1 x + 2.3 a + 8"
    @fact affToStr(aff - uaff) => "7.1 x + -2.3 a - 3"
    @fact affToStr(aff * uaff) => "(16.33 a + 39.05) x + 5.75 a + 13.75"
    @fact_throws aff / uaff
    # AffExpr--FullAffExpr
    @fact affToStr(aff + faff) => "7.1 x + (5 a + 1) x + 2 b + 5.5"
    @fact affToStr(aff - faff) => "7.1 x + (-5 a - 1) x + -2 b - 0.5"
    @fact_throws aff * faff
    @fact_throws aff / faff
end


context("Uncertain--... tests") do
    # Uncertain--Number
    @fact affToStr(a + 4.13) => "a + 4.13"
    @fact affToStr(a - 3.16) => "a - 3.16"
    @fact affToStr(a * 5.23) => "5.23 a"
    @fact affToStr(a / 2.0) => "0.5 a"
    # Uncertain--Variable
    @fact affToStr(a + x) => "x + a"
    @fact affToStr(a - x) => "-x + a"
    @fact affToStr(a * x) => "a x"
    @fact_throws affToStr(a / x)
    # Uncertain--AffExpr
    @fact affToStr(a + aff) => "7.1 x + a + 2.5"
    @fact affToStr(a - aff) => "-7.1 x + a - 2.5"
    @fact affToStr(a * aff) => "(7.1 a) x + 2.5 a"
    @fact_throws a / aff
    # Uncertain--Uncertain
    @fact affToStr(a + b) => "a + b"
    @fact affToStr(a - b) => "a - b"
    @fact_throws a * b
    @fact_throws a / b
    # Uncertain--UAffExpr (uaff = 2.3 * a + 5.5)
    @fact affToStr(b + uaff) => "b + 2.3 a + 5.5"
    @fact affToStr(b - uaff) => "b - 2.3 a - 5.5"
    @fact_throws b * uaff
    @fact_throws b / uaff
    # Uncertain--FullAffExpr (faff = (5a + 1)x + 2b + 3)
    @fact affToStr(a + faff) => "(5 a + 1) x + a + 2 b + 3"
    @fact affToStr(a - faff) => "(-5 a - 1) x + a - 2 b - 3"
    @fact_throws a * faff
    @fact_throws a / faff
end


context("UAffExpr--... tests") do
    # UAffExpr--Number
    @fact affToStr(uaff + 4.0) => "2.3 a + 9.5"
    @fact affToStr(uaff - 3.0) => "2.3 a + 2.5"
    @fact affToStr(uaff * 2.0) => "4.6 a + 11"
    @fact affToStr(uaff / 2.0) => "1.15 a + 2.75"
    # UAffExpr--Variable
    @fact affToStr(uaff + x) => "x + 2.3 a + 5.5"
    @fact affToStr(uaff - x) => "-x + 2.3 a + 5.5"
    @fact affToStr(uaff * x) => "(2.3 a + 5.5) x"
    @fact_throws uaff / x
    # UAffExpr--AffExpr (aff = 7.1 x + 2.5)
    @fact affToStr(uaff + aff) => "7.1 x + 2.3 a + 8"
    @fact affToStr(uaff - aff) => "-7.1 x + 2.3 a + 3"
    @fact affToStr(uaff * aff) => "(16.33 a + 39.05) x + 5.75 a + 13.75"
    @fact_throws uaff / aff
    # UAffExpr--Uncertain
    @fact affToStr(uaff + b) => "b + 2.3 a + 5.5"
    @fact affToStr(uaff - b) => "-b + 2.3 a + 5.5"
    @fact_throws uaff * b
    @fact_throws uaff / b
    # UAffExpr--UAffExpr (uaff2 = 3.4 b + 1.1)
    @fact affToStr(uaff + uaff2) => "2.3 a + 3.4 b + 6.6"
    @fact affToStr(uaff - uaff2) => "2.3 a - 3.4 b + 4.4"
    @fact_throws uaff * uaff2
    @fact_throws uaff / uaff2
    # UAffExpr--FullAffExpr (faff = (5a + 1)x + 2b + 3)
    @fact affToStr(uaff + faff) => "(5 a + 1) x + 2.3 a + 2 b + 8.5"
    @fact affToStr(uaff - faff) => "(-5 a - 1) x + 2.3 a - 2 b + 2.5"
    @fact_throws uaff * faff
    @fact_throws uaff / faff
end


context("FullAffExpr--... tests") do
    # faff = (5a + 1)x + 2b + 3)
    # FullAffExpr--Number
    @fact affToStr(faff + 4.0) => "(5 a + 1) x + 2 b + 7"
    @fact affToStr(faff - 2.0) => "(5 a + 1) x + 2 b + 1"
    @fact affToStr(faff * 2.0) => "(10 a + 2) x + 4 b + 6"
    @fact affToStr(faff / 2.0) => "(2.5 a + 0.5) x + b + 1.5"
    @fact conToStr(faff == 1.0) => "(5 a + 1) x + 2 b $eq -2"
    # FullAffExpr--Variable
    @fact affToStr(faff + y) => "(5 a + 1) x + y + 2 b + 3"
    @fact affToStr(faff - y) => "(5 a + 1) x - y + 2 b + 3"
    @fact_throws faff * y
    @fact_throws faff / y
    # FullAffExpr--AffExpr (aff2 = 1.2y + 1.2)
    @fact affToStr(faff + aff2) => "1.2 y + (5 a + 1) x + 2 b + 4.2"
    @fact affToStr(faff - aff2) => "(5 a + 1) x - 1.2 y + 2 b + 1.8"
    @fact_throws faff * aff2
    @fact_throws faff / aff2
    # FullAffExpr--Uncertain
    @fact affToStr(faff + a) => "(5 a + 1) x + a + 2 b + 3"
    @fact affToStr(faff - a) => "(5 a + 1) x + -a + 2 b + 3"
    @fact_throws faff * a
    @fact_throws faff / a
    # FullAffExpr--UAffExpr (uaff = 2.3 * a + 5.5)
    @fact affToStr(faff + uaff) => "(5 a + 1) x + 2.3 a + 2 b + 8.5"
    @fact affToStr(faff - uaff) => "(5 a + 1) x + 2 b - 2.3 a - 2.5"
    @fact_throws faff * uaff
    @fact_throws faff / uaff
    # FullAffExpr--FullAffExpr
    @fact affToStr(faff + faff) => "(5 a + 1) x + (5 a + 1) x + 4 b + 6"
    @fact affToStr(faff - faff) => "(5 a + 1) x + (-5 a - 1) x"
    @fact_throws faff * faff
    @fact_throws faff / faff
end

end


facts("[operators] Higher level operations") do

m = RobustModel()
@defVar(m, 0 <= x[1:3] <= 1)
@defUnc(m, 2 <= a <= 3)
@defUnc(m, 5 <= b <= 6)
@defUnc(m, u[1:3])
@defUnc(m, v[4:6])
@defUnc(m, U[1:3,1:3])
@defUnc(m, w[[:foo,:bar]])

context("sum()") do
    @fact affToStr(sum(u)) => "u[1] + u[2] + u[3]"
    @fact affToStr(sum(U)) => "U[1,1] + U[2,1] + U[3,1] + U[1,2] + U[2,2] + U[3,2] + U[1,3] + U[2,3] + U[3,3]"
    @fact affToStr(sum(w)) => anyof("w[foo] + w[bar]", "w[bar] + w[foo]")
    @fact affToStr(sum([2.0*a, 4.0*b])) => "2 a + 4 b"
    @fact affToStr(sum([x[1] + 2.0*a, x[2] + 4.0*b])) => "x[1] + x[2] + 2 a + 4 b"
end

context("dot()") do
    c = [3.5, 4.0, 2.0]
    A = [3.0 4.0 5.0;
         1.5 2.5 3.3;
         5.5 6.2 1.2]
    # DOT
    # Vector{Float64} :: JuMPArray{Uncertain}
    @fact affToStr(dot(c,u)) => "3.5 u[1] + 4 u[2] + 2 u[3]"
    @fact affToStr(dot(u,c)) => "3.5 u[1] + 4 u[2] + 2 u[3]"
    # Vector{Float64} :: JuMPArray{Uncertain} (different indices)
    @fact affToStr(dot(c,v)) => "3.5 v[4] + 4 v[5] + 2 v[6]"
    @fact affToStr(dot(v,c)) => "3.5 v[4] + 4 v[5] + 2 v[6]"
    # Array{Float64,2} (2D) :: JuMPArray{Uncertain} (2D)
    @fact affToStr(dot(A,U)) => "3 U[1,1] + 1.5 U[2,1] + 5.5 U[3,1] + 4 U[1,2] + 2.5 U[2,2] + 6.2 U[3,2] + 5 U[1,3] + 3.3 U[2,3] + 1.2 U[3,3]"

    # JuMPArray{Variable} :: JuMPArray{Uncertain}
    @fact affToStr(dot(x,u)) => "u[1] x[1] + u[2] x[2] + u[3] x[3]"

    # Array{Float64,2} (2D) :: JuMPDict{Uncertain} (1D)
    @fact_throws dot(A, u)
end

end