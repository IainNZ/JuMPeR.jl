using JuMPeR
using Base.Test

lastc(rm) = conToStr(JuMPeR.getRobust(rm).uncertainconstr[end])
lastuc(rm) = conToStr(JuMPeR.getRobust(rm).uncertaintyset[end])

rm = RobustModel()
@defVar(rm, x)
@defVar(rm, y[1:5])
@defUnc(rm, u)
@defUnc(rm, v[1:5])

@addConstraint(rm, u*x <= 5)
@test lastc(rm) == "u x <= 5"

@addConstraint(rm, sum{v[i]*y[i], i=1:5; i!=3} <= 9)
@test lastc(rm) == "v[1] y[1] + v[2] y[2] + v[4] y[4] + v[5] y[5] <= 9"

@addConstraint(rm, sum{v[i], i=1:5} == 1)
@test lastuc(rm) == "v[5] + v[4] + v[3] + v[2] + v[1] == 1"

@addConstraint(rm, u == sum{i*v[i], i=1:3})
@test lastuc(rm) == "-3 v[3] - 2 v[2] - v[1] + u == 0"

@addConstraint(rm, sum{i*(u+v[i])*(y[i]+x), i=1:2:5} <= 0)
@test lastc(rm) == "(u + v[1]) y[1] + (u + v[1]) x + (3 u + 3 v[3]) y[3] + (3 u + 3 v[3]) x + (5 u + 5 v[5]) y[5] + (5 u + 5 v[5]) x <= 0"