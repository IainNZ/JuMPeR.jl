using JuMPeR
using FactCheck, Base.Test

facts("[operators] Robust operator tests") do

m = RobustModel()
@defVar(m, w)
@defVar(m, x)
@defVar(m, y)
@defVar(m, z)

aff = 7.1 * x + 2.5
@fact affToStr(aff) => "7.1 x + 2.5"

aff2 = 1.2 * y + 1.2
@fact affToStr(aff2) => "1.2 y + 1.2"

a = Uncertain(m, 2.0, 3.0, :Cont, "a")
b = Uncertain(m, 5.0, 6.0, :Cont, "b")

uaff = 2.3 * a + 5.5
@fact affToStr(uaff) => "2.3 a + 5.5"
uaff2 = 3.4 * b + 1.1
@fact affToStr(uaff2) => "3.4 b + 1.1"

faff = FullAffExpr([x],[UAffExpr([a],[5.],1.)],UAffExpr([b],[2.],3.))
@fact affToStr(faff) => "(5 a + 1) x + 2 b + 3"


context("Number--... tests") do
    # Number--Uncertain
    @fact  affToStr(4.13 + a) => "a + 4.13"
    @fact  affToStr(3.16 - a) => "-a + 3.16"
    @fact  affToStr(5.23 * a) => "5.23 a"
    @fact_throws 2.94 / a
    # Number--UAffExpr
    @fact  affToStr(2.3 + uaff) => "2.3 a + 7.8"
    @fact  affToStr(1.5 - uaff) => "-2.3 a - 4"
    @fact  affToStr(2.0 * uaff) => "4.6 a + 11"
    @fact_throws 2.94 / uaff
    # Number--FullAffExpr
    @fact  affToStr(2.3 + faff) => "(5 a + 1) x + 2 b + 5.3"
    @fact  affToStr(1.0 - faff) => "(-5 a - 1) x + -2 b - 2"
    @fact  affToStr(2.0 * faff) => "(10 a + 2) x + 4 b + 6"
    @fact_throws 2.94 / faff
end


context("Variable--... tests") do
    # Variable--Uncertain
    @fact  affToStr(x + a) => "x + a"
    @fact  affToStr(x - a) => "x + -a"
    @fact  affToStr(x * a) => "a x"
    @fact_throws affToStr(x / a)
    # Variable--UAffExpr
    @fact  affToStr(x + uaff) => "x + 2.3 a + 5.5"
    @fact  affToStr(x - uaff) => "x + -2.3 a - 5.5"
    @fact  affToStr(x * uaff) => "(2.3 a + 5.5) x"
    @fact_throws affToStr(x / uaff)
    # Variable--FullAffExpr
    @fact  affToStr(x + faff) => "(5 a + 1) x + x + 2 b + 3"
    @fact  affToStr(x - faff) => "(-5 a - 1) x + x + -2 b - 3"
    @fact_throws x * faff
    @fact_throws x / faff
end


context("AffExpr--... tests") do
    # AffExpr--Uncertain
    @fact  affToStr(aff + a) => "7.1 x + a + 2.5"
    @fact  affToStr(aff - a) => "7.1 x + -a + 2.5"
    @fact  affToStr(aff * a) => "(7.1 a) x + 2.5 a"
    @fact_throws aff / a
    # AffExpr--UAffExpr
    @fact  affToStr(aff + uaff) => "7.1 x + 2.3 a + 8"
    @fact  affToStr(aff - uaff) => "7.1 x + -2.3 a - 3"
    @fact  affToStr(aff * uaff) => "(16.33 a + 39.05) x + 5.75 a + 13.75"
    @fact_throws aff / uaff
    # AffExpr--FullAffExpr
    @fact  affToStr(aff + faff) => "7.1 x + (5 a + 1) x + 2 b + 5.5"
    @fact  affToStr(aff - faff) => "7.1 x + (-5 a - 1) x + -2 b - 0.5"
    @fact_throws aff * faff
    @fact_throws aff / faff
end


context("Uncertain--... tests") do
    # Uncertain--Number
    @fact  affToStr(a + 4.13) => "a + 4.13"
    @fact  affToStr(a - 3.16) => "a - 3.16"
    @fact  affToStr(a * 5.23) => "5.23 a"
    @fact  affToStr(a / 2.0) => "0.5 a"
    # Uncertain--Variable
    @fact  affToStr(a + x) => "x + a"
    @fact  affToStr(a - x) => "-x + a"
    @fact  affToStr(a * x) => "a x"
    @fact_throws affToStr(a / x)
    # Uncertain--AffExpr
    @fact  affToStr(a + aff) => "7.1 x + a + 2.5"
    @fact  affToStr(a - aff) => "-7.1 x + a - 2.5"
    @fact  affToStr(a * aff) => "(7.1 a) x + 2.5 a"
    @fact_throws a / aff
    # Uncertain--Uncertain
    @fact  affToStr(a + b) => "a + b"
    @fact  affToStr(a - b) => "a - b"
    @fact_throws a * b
    @fact_throws a / b
    # Uncertain--UAffExpr (uaff = 2.3 * a + 5.5)
    @fact  affToStr(b + uaff) => "b + 2.3 a + 5.5"
    @fact  affToStr(b - uaff) => "b - 2.3 a - 5.5"
    @fact_throws b * uaff
    @fact_throws b / uaff
    # Uncertain--FullAffExpr (faff = (5a + 1)x + 2b + 3)
    @fact  affToStr(a + faff) => "(5 a + 1) x + a + 2 b + 3"
    @fact  affToStr(a - faff) => "(-5 a - 1) x + a - 2 b - 3"
    @fact_throws a * faff
    @fact_throws b * faff
end


context("UAffExpr--... tests") do
    # UAffExpr--Number
    @fact  affToStr(uaff + 4.0) => "2.3 a + 9.5"
    @fact  affToStr(uaff - 3.0) => "2.3 a + 2.5"
    @fact  affToStr(uaff * 2.0) => "4.6 a + 11"
    @fact  affToStr(uaff / 2.0) => "1.15 a + 2.75"
    # UAffExpr--Variable
    @fact  affToStr(uaff + x) => "x + 2.3 a + 5.5"
    @fact  affToStr(uaff - x) => "-x + 2.3 a + 5.5"
    @fact  affToStr(uaff * x) => "(2.3 a + 5.5) x"
    @fact_throws uaff / x
    # UAffExpr--AffExpr (aff = 7.1 x + 2.5)
    @fact  affToStr(uaff + aff) => "7.1 x + 2.3 a + 8"
    @fact  affToStr(uaff - aff) => "-7.1 x + 2.3 a + 3"
    @fact  affToStr(uaff * aff) => "(16.33 a + 39.05) x + 5.75 a + 13.75"
    @fact_throws uaff / aff
    # UAffExpr--Uncertain
    @fact  affToStr(uaff + b) => "b + 2.3 a + 5.5"
    @fact  affToStr(uaff - b) => "-b + 2.3 a + 5.5"
    @fact_throws uaff * b
    @fact_throws uaff / b
    # UAffExpr--UAffExpr (uaff2 = 3.4 b + 1.1)
    @fact  affToStr(uaff + uaff2) => "2.3 a + 3.4 b + 6.6"
    @fact  affToStr(uaff - uaff2) => "2.3 a - 3.4 b + 4.4"
    @fact_throws uaff * uaff2
    @fact_throws uaff / uaff2
    # UAffExpr--FullAffExpr (faff = (5a + 1)x + 2b + 3)
    @fact  affToStr(uaff + faff) => "(5 a + 1) x + 2.3 a + 2 b + 8.5"
    @fact  affToStr(uaff - faff) => "(-5 a - 1) x + 2.3 a - 2 b + 2.5"
    @fact_throws uaff * faff
    @fact_throws uaff / faff
end


context("FullAffExpr--... tests") do
    # faff = (5a + 1)x + 2b + 3)
    # FullAffExpr--Number
    @fact  affToStr(faff + 4.0) => "(5 a + 1) x + 2 b + 7"
    @fact  affToStr(faff - 2.0) => "(5 a + 1) x + 2 b + 1"
    @fact  affToStr(faff * 2.0) => "(10 a + 2) x + 4 b + 6"
    @fact  affToStr(faff / 2.0) => "(2.5 a + 0.5) x + b + 1.5"
    # FullAffExpr--Variable
    @fact  affToStr(faff + y) => "(5 a + 1) x + y + 2 b + 3"
    @fact  affToStr(faff - y) => "(5 a + 1) x - y + 2 b + 3"
    @fact_throws faff * y
    @fact_throws faff / y
    # FullAffExpr--AffExpr (aff2 = 1.2y + 1.2)
    @fact  affToStr(faff + aff2) => "1.2 y + (5 a + 1) x + 2 b + 4.2"
    @fact  affToStr(faff - aff2) => "(5 a + 1) x - 1.2 y + 2 b + 1.8"
    @fact_throws faff * aff2
    @fact_throws faff / aff2
    # FullAffExpr--Uncertain
    @fact  affToStr(faff + a) => "(5 a + 1) x + a + 2 b + 3"
    @fact  affToStr(faff - a) => "(5 a + 1) x + -a + 2 b + 3"
    @fact_throws faff * a
    @fact_throws faff / a
    # FullAffExpr--UAffExpr (uaff = 2.3 * a + 5.5)
    @fact  affToStr(faff + uaff) => "(5 a + 1) x + 2.3 a + 2 b + 8.5"
    @fact  affToStr(faff - uaff) => "(5 a + 1) x + 2 b - 2.3 a - 2.5"
    @fact_throws faff * uaff
    @fact_throws faff / uaff
    # FullAffExpr--FullAffExpr
    @fact  affToStr(faff + faff) => "(5 a + 1) x + (5 a + 1) x + 4 b + 6"
    @fact  affToStr(faff - faff) => "(5 a + 1) x + (-5 a - 1) x"
    @fact_throws faff * faff
    @fact_throws faff / faff
end

end

# Higher level
m2 = RobustModel()
@defVar(m2, w)
@defVar(m2, x)
@defVar(m2, y)
@defVar(m2, z)
a = Uncertain(m2, 2.0, 3.0, :Cont, "a")
b = Uncertain(m2, 5.0, 6.0, :Cont, "b")
@defUnc(m2, udict[1:3])
setName(udict[1], "a")
setName(udict[2], "b")
setName(udict[3], "c")
@defUnc(m2, vdict[4:6])
setName(vdict[4], "d")
setName(vdict[5], "e")
setName(vdict[6], "f")
@defUnc(m2, matdict[1:3,1:3])
for i = 1:3
    for j = 1:3
        setName(matdict[i,j], "U$(i)$(j)")
    end
end
@defVar(m2, 0 <= x[1:3] <= 1)

nums = [3.5, 4.0, 2.0]
A = [3.0 4.0 5.0;
     1.5 2.5 3.3;
     5.5 6.2 1.2]

# SUM
@fact affToStr(sum(udict)) => "a + b + c"
@fact affToStr(sum(matdict)) => "U11 + U21 + U31 + U12 + U22 + U32 + U13 + U23 + U33"
@fact affToStr(sum([2.0*a, 4.0*b])) => "2 a + 4 b"
@fact affToStr(sum([x[1] + 2.0*a, x[2] + 4.0*b])) => "x[1] + x[2] + 2 a + 4 b"

# DOT
# JuMPDict{Variable} :: JuMPDict{Uncertain}
#@fact affToStr(dot(x, udict)) => "a x[1] + b x[2] + c x[3]"  # TODO
# Vector{Float64} :: JuMPDict{Uncertain}
@fact affToStr(dot(nums,udict)) => "3.5 a + 4 b + 2 c"
@fact affToStr(dot(udict,nums)) => "3.5 a + 4 b + 2 c"
# Vector{Float64} :: JuMPDict{Uncertain} (different indices)
@fact affToStr(dot(nums,vdict)) => "3.5 d + 4 e + 2 f"
@fact affToStr(dot(vdict,nums)) => "3.5 d + 4 e + 2 f"
# Array{Float64,2} (1D) :: JuMPDict{Uncertain} 
@fact affToStr(dot(vec(A[1,:]), udict)) => "3 a + 4 b + 5 c"
@fact affToStr(dot(vec(A[1,:]), vdict)) => "3 d + 4 e + 5 f"
# Array{Float64,2} (2D) :: JuMPDict{Uncertain} (2D)
@fact affToStr(dot(A, matdict)) => "3 U11 + 1.5 U21 + 5.5 U31 + 4 U12 + 2.5 U22 + 6.2 U32 + 5 U13 + 3.3 U23 + 1.2 U33"
# Vector{UAffExpr} :: JuMPDict{Variable}
@fact affToStr(dot(
    [1*udict[1],2*udict[2],3*udict[3]],
    x)) => "a x[1] + (2 b) x[2] + (3 c) x[3]"
# Array{Float64,2} (2D) :: JuMPDict{Uncertain} (1D)
@fact_throws dot(A, udict)
# Array{Float64,1} (1D) :: JuMPDict{Uncertain} (2D)
@fact_throws dot(nums, matdict)