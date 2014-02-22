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
@test affToStr(faff) == "(5.0 a + 1.0) x + 2.0 b + 3.0"

# 1. Number tests
# Number--Uncertain
@test affToStr(4.13 + a) == "1.0 a + 4.13"
@test affToStr(3.16 - a) == "-1.0 a + 3.16"
@test affToStr(5.23 * a) == "5.23 a"
@test_throws 2.94 / a
# Number--UAffExpr
@test affToStr(2.3 + uaff) == "2.3 a + 7.8"
@test affToStr(1.5 - uaff) == "-2.3 a + -4.0"
@test affToStr(2.0 * uaff) == "4.6 a + 11.0"
@test_throws 2.94 / uaff
# Number--FullAffExpr
@test affToStr(2.3 + faff) == "(5.0 a + 1.0) x + 2.0 b + 5.3"
@test affToStr(1.0 - faff) == "(-5.0 a + -1.0) x + -2.0 b + -2.0"
@test affToStr(2.0 * faff) == "(10.0 a + 2.0) x + 4.0 b + 6.0"
@test_throws 2.94 / faff

# 2. Variable test
# Variable--Uncertain
@test affToStr(x + a) == "(1.0) x + 1.0 a"
@test affToStr(x - a) == "(1.0) x + -1.0 a"
@test affToStr(x * a) == "(1.0 a) x + 0.0"
@test_throws affToStr(x / a)
# Variable--UAffExpr
@test affToStr(x + uaff) == "(1.0) x + 2.3 a + 5.5"
@test affToStr(x - uaff) == "(1.0) x + -2.3 a + -5.5"
@test affToStr(x * uaff) == "(2.3 a + 5.5) x + 0.0"
@test_throws affToStr(x / uaff)
# Variable--FullAffExpr
@test affToStr(x + faff) == "(5.0 a + 1.0) x + (1.0) x + 2.0 b + 3.0"
@test affToStr(x - faff) == "(-5.0 a + -1.0) x + (1.0) x + -2.0 b + -3.0"
@test_throws x * faff
@test_throws x / faff

# 3. AffExpr test
# AffExpr--Uncertain
@test affToStr(aff + a) == "(7.1) x + 1.0 a + 2.5"
@test affToStr(aff - a) == "(7.1) x + -1.0 a + 2.5"
@test affToStr(aff * a) == "(7.1 a) x + 2.5 a"
@test_throws aff / a
# AffExpr--UAffExpr
@test affToStr(aff + uaff) == "(7.1) x + 2.3 a + 8.0"
@test affToStr(aff - uaff) == "(7.1) x + -2.3 a + -3.0"
@test affToStr(aff * uaff) == "(16.33 a + 39.05) x + 5.75 a + 13.75"
@test_throws aff / uaff
# AffExpr--FullAffExpr
@test affToStr(aff + faff) == "(7.1) x + (5.0 a + 1.0) x + 2.0 b + 5.5"
@test affToStr(aff - faff) == "(7.1) x + (-5.0 a + -1.0) x + -2.0 b + -0.5"
@test_throws aff * faff
@test_throws aff / faff

# 6. Uncertain test
# Uncertain--Number
@test affToStr(a + 4.13) == "1.0 a + 4.13"
@test affToStr(a - 3.16) == "1.0 a + -3.16"
@test affToStr(a * 5.23) == "5.23 a"
@test affToStr(a / 2.0) == "0.5 a"
# Uncertain--Variable
@test affToStr(a + x) == "(1.0) x + 1.0 a"
@test affToStr(a - x) == "(-1.0) x + 1.0 a"
@test affToStr(a * x) == "(1.0 a) x + 0.0"
@test_throws affToStr(a / x)
# Uncertain--AffExpr
@test affToStr(a + aff) == "(7.1) x + 1.0 a + 2.5"
@test affToStr(a - aff) == "(-7.1) x + 1.0 a + -2.5"
@test affToStr(a * aff) == "(7.1 a) x + 2.5 a"
@test_throws a / aff
# Uncertain--Uncertain
@test affToStr(a + b) == "1.0 a + 1.0 b"
@test affToStr(a - b) == "1.0 a + -1.0 b"
@test_throws a * b
@test_throws a / b
# Uncertain--UAffExpr (uaff = 2.3 * a + 5.5)
@test affToStr(b + uaff) == "1.0 b + 2.3 a + 5.5"
@test affToStr(b - uaff) == "1.0 b + -2.3 a + -5.5"
@test_throws b * uaff
@test_throws b / uaff
# Uncertain--FullAffExpr (faff = (5a + 1)x + 2b + 3)
@test affToStr(a + faff) == "(5.0 a + 1.0) x + 1.0 a + 2.0 b + 3.0"
@test affToStr(a - faff) == "(-5.0 a + -1.0) x + 1.0 a + -2.0 b + -3.0"
@test_throws a * faff
@test_throws b * faff

# 7. UAffExpr test (uaff = 2.3 * a + 5.5)
# UAffExpr--Number
@test affToStr(uaff + 4.0) == "2.3 a + 9.5"
@test affToStr(uaff - 3.0) == "2.3 a + 2.5"
@test affToStr(uaff * 2.0) == "4.6 a + 11.0"
@test affToStr(uaff / 2.0) == "1.15 a + 2.75"
# UAffExpr--Variable
@test affToStr(uaff + x) == "(1.0) x + 2.3 a + 5.5"
@test affToStr(uaff - x) == "(-1.0) x + 2.3 a + 5.5"
@test affToStr(uaff * x) == "(2.3 a + 5.5) x + 0.0"
@test_throws uaff / x
# UAffExpr--AffExpr (aff = 7.1 x + 2.5)
@test affToStr(uaff + aff) == "(7.1) x + 2.3 a + 8.0"
@test affToStr(uaff - aff) == "(-7.1) x + 2.3 a + 3.0"
@test affToStr(uaff * aff) == "(16.33 a + 39.05) x + 5.75 a + 13.75"
@test_throws uaff / aff
# UAffExpr--Uncertain
@test affToStr(uaff + b) == "1.0 b + 2.3 a + 5.5"
@test affToStr(uaff - b) == "-1.0 b + 2.3 a + 5.5"
@test_throws uaff * b
@test_throws uaff / b
# UAffExpr--UAffExpr (uaff2 = 3.4 b + 1.1)
@test affToStr(uaff + uaff2) == "2.3 a + 3.4 b + 6.6"
@test affToStr(uaff - uaff2) == "2.3 a + -3.4 b + 4.4"
@test_throws uaff * uaff2
@test_throws uaff / uaff2
# UAffExpr--FullAffExpr (faff = (5a + 1)x + 2b + 3)
@test affToStr(uaff + faff) == "(5.0 a + 1.0) x + 2.3 a + 2.0 b + 8.5"
@test affToStr(uaff - faff) == "(-5.0 a + -1.0) x + 2.3 a + -2.0 b + 2.5"
@test_throws uaff * faff
@test_throws uaff / faff

# 8. FullAffExpr test (faff = (5a + 1)x + 2b + 3)
# FullAffExpr--Number
@test affToStr(faff + 4.0) == "(5.0 a + 1.0) x + 2.0 b + 7.0"
@test affToStr(faff - 2.0) == "(5.0 a + 1.0) x + 2.0 b + 1.0"
@test affToStr(faff * 2.0) == "(10.0 a + 2.0) x + 4.0 b + 6.0"
@test affToStr(faff / 2.0) == "(2.5 a + 0.5) x + 1.0 b + 1.5"
# FullAffExpr--Variable
@test affToStr(faff + y) == "(5.0 a + 1.0) x + (1.0) y + 2.0 b + 3.0"
@test affToStr(faff - y) == "(5.0 a + 1.0) x + (-1.0) y + 2.0 b + 3.0"
@test_throws faff * y
@test_throws faff / y
# FullAffExpr--AffExpr (aff2 = 1.2y + 1.2)
@test affToStr(faff + aff2) == "(1.2) y + (5.0 a + 1.0) x + 2.0 b + 4.2"
@test affToStr(faff - aff2) == "(5.0 a + 1.0) x + (-1.2) y + 2.0 b + 1.8"
@test_throws faff * aff2
@test_throws faff / aff2
# FullAffExpr--Uncertain
@test affToStr(faff + a) == "(5.0 a + 1.0) x + 1.0 a + 2.0 b + 3.0"
@test affToStr(faff - a) == "(5.0 a + 1.0) x + -1.0 a + 2.0 b + 3.0"
@test_throws faff * a
@test_throws faff / a
# FullAffExpr--UAffExpr (uaff = 2.3 * a + 5.5)
@test affToStr(faff + uaff) == "(5.0 a + 1.0) x + 2.3 a + 2.0 b + 8.5"
@test affToStr(faff - uaff) == "(5.0 a + 1.0) x + 2.0 b + -2.3 a + -2.5"
@test_throws faff * uaff
@test_throws faff / uaff
# FullAffExpr--FullAffExpr
@test affToStr(faff + faff) == "(5.0 a + 1.0) x + (5.0 a + 1.0) x + 4.0 b + 6.0"
@test affToStr(faff - faff) == "(5.0 a + 1.0) x + (-5.0 a + -1.0) x + 1.0e-50 b"
@test_throws faff * faff
@test_throws faff / faff
