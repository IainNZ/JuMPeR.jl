#############################################################################
# JuMPeR
# Julia for Mathematical Programming - extension for Robust Optimization
# See http://github.com/IainNZ/JuMPeR.jl
#############################################################################
# test/macro.jl
# Testing JuMP-defined macros still work with Uncertains etc.
#############################################################################

using JuMPeR
using FactCheck

lastuc(rm) = conToStr(JuMPeR.getRobust(rm).uncertainconstr[end])
lastus(rm) = conToStr(JuMPeR.getRobust(rm).uncertaintyset[end])
le, ge = JuMP.repl_leq, JuMP.repl_geq


facts("[macro] Uncertainty set constraints") do

    rm = RobustModel()
    @defUnc(rm, u)
    @defUnc(rm, v[1:3])

    @addConstraint(rm, u <= 5)
    @fact lastus(rm) => "u $le 5"
    @addConstraint(rm, u - 5 <= 0)
    @fact lastus(rm) => "u $le 5"
    @addConstraint(rm, 5 <= u)
    @fact lastus(rm) => "-u $le -5"
    @addConstraint(rm, 0 <= u - 5)
    @fact lastus(rm) => "-u $le -5"

    @addConstraint(rm, u == sum{i*v[i], i=1:3})
    @fact lastus(rm) => "-3 v[3] - 2 v[2] - v[1] + u = 0"
    @addConstraint(rm, sum{i*v[i], i=1:3} + u >= 10)
    @fact lastus(rm) => "u + 3 v[3] + 2 v[2] + v[1] $ge 10"
end

facts("[macro] Uncertain constraints") do

rm = RobustModel()
@defVar(rm, x)
@defVar(rm, y[1:3])
@defVar(rm, z)

@defUnc(rm, u)
@defUnc(rm, v[1:3])
@defUnc(rm, w)

context("Simple u*x") do
    @addConstraint(rm, u*x <= 5)
    @fact lastuc(rm) => "u x $le 5"

    @addConstraint(rm, x*u <= 5)
    @fact lastuc(rm) => "u x $le 5"

    @addConstraint(rm, 5 <= u*x)
    @fact lastuc(rm) => "-u x $le -5"

    @addConstraint(rm, 5 <= x*u)
    @fact lastuc(rm) => "-u x $le -5"
end

context("Adding fullaffs") do
    @addConstraint(rm, (u+w)*x + 2 + w*x <= u*z + 3)
    @fact lastuc(rm) => "(u + w) x + w x + -u z $le 1"
end

end



#=


@addConstraint(rm, sum{v[i]*y[i], i=1:5; i!=3} <= 9)
@test lastc(rm) == "v[1] y[1] + v[2] y[2] + v[4] y[4] + v[5] y[5] <= 9"

@addConstraint(rm, sum{v[i], i=1:5} == 1)
@test lastuc(rm) == "v[5] + v[4] + v[3] + v[2] + v[1] == 1"


@addConstraint(rm, sum{i*(u+v[i])*(y[i]+x), i=1:2:5} <= 0)
@test lastc(rm) == "(u + v[1]) y[1] + (u + v[1]) x + (3 u + 3 v[3]) y[3] + (3 u + 3 v[3]) x + (5 u + 5 v[5]) y[5] + (5 u + 5 v[5]) x <= 0"

end
=#