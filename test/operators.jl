using JuMPeR
using Base.Test

m = RobustModel()
@defVar(m, w)
@defVar(m, x)
@defVar(m, y)
@defVar(m, z)

aff = 7.1 * x + 2.5
@test affToStr(aff) == "7.1 x + 2.5"

aff2 = 1.2 * y + 1.2
@test affToStr(aff2) == "1.2 y + 1.2"

a = Uncertain(m, 2.0, 3.0, "a")
b = Uncertain(m, 5.0, 6.0, "b")

uaff = 2.3 * a + 5.5
@test affToStr(uaff) == "2.3 a + 5.5"
uaff2 = 3.4 * b + 1.1
@test affToStr(uaff2) == "3.4 b + 1.1"

faff = FullAffExpr([x],[UAffExpr([a],[5.],1.)],UAffExpr([b],[2.],3.))
@test affToStr(faff) == "(5 a + 1) x + 2 b + 3"

# 1. Number tests
# Number--Uncertain
@test affToStr(4.13 + a) == "a + 4.13"
@test affToStr(3.16 - a) == "-a + 3.16"
@test affToStr(5.23 * a) == "5.23 a"
@test_throws    2.94 / a
# Number--UAffExpr
@test affToStr(2.3 + uaff) == "2.3 a + 7.8"
@test affToStr(1.5 - uaff) == "-2.3 a - 4"
@test affToStr(2.0 * uaff) == "4.6 a + 11"
@test_throws    2.94 / uaff
# Number--FullAffExpr
@test affToStr(2.3 + faff) == "(5 a + 1) x + 2 b + 5.3"
@test affToStr(1.0 - faff) == "(-5 a - 1) x + -2 b - 2"
@test affToStr(2.0 * faff) == "(10 a + 2) x + 4 b + 6"
@test_throws    2.94 / faff

# 2. Variable test
# Variable--Uncertain
@test affToStr(x + a) == "x + a"
@test affToStr(x - a) == "x + -a"
@test affToStr(x * a) == "a x"
@test_throws    affToStr(x / a)
# Variable--UAffExpr
@test affToStr(x + uaff) == "x + 2.3 a + 5.5"
@test affToStr(x - uaff) == "x + -2.3 a - 5.5"
@test affToStr(x * uaff) == "(2.3 a + 5.5) x"
@test_throws    affToStr(x / uaff)
# Variable--FullAffExpr
@test affToStr(x + faff) == "(5 a + 1) x + x + 2 b + 3"
@test affToStr(x - faff) == "(-5 a - 1) x + x + -2 b - 3"
@test_throws    x * faff
@test_throws    x / faff

# 3. AffExpr test
# AffExpr--Uncertain
@test affToStr(aff + a) == "7.1 x + a + 2.5"
@test affToStr(aff - a) == "7.1 x + -a + 2.5"
@test affToStr(aff * a) == "(7.1 a) x + 2.5 a"
@test_throws    aff / a
# AffExpr--UAffExpr
@test affToStr(aff + uaff) == "7.1 x + 2.3 a + 8"
@test affToStr(aff - uaff) == "7.1 x + -2.3 a - 3"
@test affToStr(aff * uaff) == "(16.33 a + 39.05) x + 5.75 a + 13.75"
@test_throws    aff / uaff
# AffExpr--FullAffExpr
@test affToStr(aff + faff) == "7.1 x + (5 a + 1) x + 2 b + 5.5"
@test affToStr(aff - faff) == "7.1 x + (-5 a - 1) x + -2 b - 0.5"
@test_throws    aff * faff
@test_throws    aff / faff

# 6. Uncertain test
# Uncertain--Number
@test affToStr(a + 4.13) == "a + 4.13"
@test affToStr(a - 3.16) == "a - 3.16"
@test affToStr(a * 5.23) == "5.23 a"
@test affToStr(a / 2.0) == "0.5 a"
# Uncertain--Variable
@test affToStr(a + x) == "x + a"
@test affToStr(a - x) == "-x + a"
@test affToStr(a * x) == "a x"
@test_throws    affToStr(a / x)
# Uncertain--AffExpr
@test affToStr(a + aff) == "7.1 x + a + 2.5"
@test affToStr(a - aff) == "-7.1 x + a - 2.5"
@test affToStr(a * aff) == "(7.1 a) x + 2.5 a"
@test_throws    a / aff
# Uncertain--Uncertain
@test affToStr(a + b) == "a + b"
@test affToStr(a - b) == "a - b"
@test_throws    a * b
@test_throws    a / b
# Uncertain--UAffExpr (uaff = 2.3 * a + 5.5)
@test affToStr(b + uaff) == "b + 2.3 a + 5.5"
@test affToStr(b - uaff) == "b - 2.3 a - 5.5"
@test_throws    b * uaff
@test_throws    b / uaff
# Uncertain--FullAffExpr (faff = (5a + 1)x + 2b + 3)
@test affToStr(a + faff) == "(5 a + 1) x + a + 2 b + 3"
@test affToStr(a - faff) == "(-5 a - 1) x + a - 2 b - 3"
@test_throws    a * faff
@test_throws    b * faff

# 7. UAffExpr test (uaff = 2.3 * a + 5.5)
# UAffExpr--Number
@test affToStr(uaff + 4.0) == "2.3 a + 9.5"
@test affToStr(uaff - 3.0) == "2.3 a + 2.5"
@test affToStr(uaff * 2.0) == "4.6 a + 11"
@test affToStr(uaff / 2.0) == "1.15 a + 2.75"
# UAffExpr--Variable
@test affToStr(uaff + x) == "x + 2.3 a + 5.5"
@test affToStr(uaff - x) == "-x + 2.3 a + 5.5"
@test affToStr(uaff * x) == "(2.3 a + 5.5) x"
@test_throws    uaff / x
# UAffExpr--AffExpr (aff = 7.1 x + 2.5)
@test affToStr(uaff + aff) == "7.1 x + 2.3 a + 8"
@test affToStr(uaff - aff) == "-7.1 x + 2.3 a + 3"
@test affToStr(uaff * aff) == "(16.33 a + 39.05) x + 5.75 a + 13.75"
@test_throws    uaff / aff
# UAffExpr--Uncertain
@test affToStr(uaff + b) == "b + 2.3 a + 5.5"
@test affToStr(uaff - b) == "-b + 2.3 a + 5.5"
@test_throws    uaff * b
@test_throws    uaff / b
# UAffExpr--UAffExpr (uaff2 = 3.4 b + 1.1)
@test affToStr(uaff + uaff2) == "2.3 a + 3.4 b + 6.6"
@test affToStr(uaff - uaff2) == "2.3 a - 3.4 b + 4.4"
@test_throws    uaff * uaff2
@test_throws    uaff / uaff2
# UAffExpr--FullAffExpr (faff = (5a + 1)x + 2b + 3)
@test affToStr(uaff + faff) == "(5 a + 1) x + 2.3 a + 2 b + 8.5"
@test affToStr(uaff - faff) == "(-5 a - 1) x + 2.3 a - 2 b + 2.5"
@test_throws    uaff * faff
@test_throws    uaff / faff

# 8. FullAffExpr test (faff = (5a + 1)x + 2b + 3)
# FullAffExpr--Number
@test affToStr(faff + 4.0) == "(5 a + 1) x + 2 b + 7"
@test affToStr(faff - 2.0) == "(5 a + 1) x + 2 b + 1"
@test affToStr(faff * 2.0) == "(10 a + 2) x + 4 b + 6"
@test affToStr(faff / 2.0) == "(2.5 a + 0.5) x + b + 1.5"
# FullAffExpr--Variable
@test affToStr(faff + y) == "(5 a + 1) x + y + 2 b + 3"
@test affToStr(faff - y) == "(5 a + 1) x - y + 2 b + 3"
@test_throws    faff * y
@test_throws    faff / y
# FullAffExpr--AffExpr (aff2 = 1.2y + 1.2)
@test affToStr(faff + aff2) == "1.2 y + (5 a + 1) x + 2 b + 4.2"
@test affToStr(faff - aff2) == "(5 a + 1) x - 1.2 y + 2 b + 1.8"
@test_throws    faff * aff2
@test_throws    faff / aff2
# FullAffExpr--Uncertain
@test affToStr(faff + a) == "(5 a + 1) x + a + 2 b + 3"
@test affToStr(faff - a) == "(5 a + 1) x + -a + 2 b + 3"
@test_throws    faff * a
@test_throws    faff / a
# FullAffExpr--UAffExpr (uaff = 2.3 * a + 5.5)
@test affToStr(faff + uaff) == "(5 a + 1) x + 2.3 a + 2 b + 8.5"
@test affToStr(faff - uaff) == "(5 a + 1) x + 2 b - 2.3 a - 2.5"
@test_throws    faff * uaff
@test_throws    faff / uaff
# FullAffExpr--FullAffExpr
@test affToStr(faff + faff) == "(5 a + 1) x + (5 a + 1) x + 4 b + 6"
@test affToStr(faff - faff) == "(5 a + 1) x + (-5 a - 1) x"
@test_throws    faff * faff
@test_throws    faff / faff

# Expression operations
push!(faff, 1, z)

# Higher level
m2 = RobustModel()
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
@test affToStr(sum(udict)) == "a + b + c"
@test affToStr(sum(matdict)) == "U11 + U21 + U31 + U12 + U22 + U32 + U13 + U23 + U33"
@test affToStr(sum([2.0*a, 4.0*b])) == "2 a + 4 b"
@test affToStr(sum([x[1] + 2.0*a, x[2] + 4.0*b])) == "x[1] + x[2] + 2 a + 4 b"

# DOT
# JuMPDict{Variable} :: JuMPDict{Uncertain}
@test affToStr(dot(x, udict)) == "a x[1] + b x[2] + c x[3]"
# Vector{Float64} :: JuMPDict{Uncertain}
@test affToStr(dot(nums,udict)) == "3.5 a + 4 b + 2 c"
@test affToStr(dot(udict,nums)) == "3.5 a + 4 b + 2 c"
# Vector{Float64} :: JuMPDict{Uncertain} (different indices)
@test affToStr(dot(nums,vdict)) == "3.5 d + 4 e + 2 f"
@test affToStr(dot(vdict,nums)) == "3.5 d + 4 e + 2 f"
# Array{Float64,2} (1D) :: JuMPDict{Uncertain} 
@test affToStr(dot(vec(A[1,:]), udict)) == "3 a + 4 b + 5 c"
@test affToStr(dot(vec(A[1,:]), vdict)) == "3 d + 4 e + 5 f"
# Array{Float64,2} (2D) :: JuMPDict{Uncertain} (2D)
@test affToStr(dot(A, matdict)) == "3 U11 + 1.5 U21 + 5.5 U31 + 4 U12 + 2.5 U22 + 6.2 U32 + 5 U13 + 3.3 U23 + 1.2 U33"
# Vector{UAffExpr} :: JuMPDict{Variable}
@test affToStr(dot(
    [1*udict[1],2*udict[2],3*udict[3]],
    x)) == "a x[1] + (2 b) x[2] + (3 c) x[3]"
# Array{Float64,2} (2D) :: JuMPDict{Uncertain} (1D)
@test_throws  dot(A, udict)
# Array{Float64,1} (1D) :: JuMPDict{Uncertain} (2D)
@test_throws  dot(nums, matdict)