#-----------------------------------------------------------------------
# JuMPeR  --  JuMP Extension for Robust Optimization
# http://github.com/IainNZ/JuMPeR.jl
#-----------------------------------------------------------------------
# Copyright (c) 2015: Iain Dunning
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
#-----------------------------------------------------------------------
# test/macro.jl
# Testing JuMP-defined macros still work with Uncertains etc.
#-----------------------------------------------------------------------

using JuMP, JuMPeR
using FactCheck

lastuc(rm) = conToStr(JuMPeR.getRobust(rm).uncertainconstr[end])
lastus(rm) = conToStr(JuMPeR.getRobust(rm).uncertaintyset[end])
le, eq, ge = JuMP.repl[:leq], JuMP.repl[:eq], JuMP.repl[:geq]


facts("[macro] Uncertainty set constraints") do

    rm = RobustModel()
    @defUnc(rm, u)
    @defUnc(rm, v[1:3])

    @addConstraint(rm, u <= 5)
    @fact lastus(rm) --> "u $le 5"
    @addConstraint(rm, u - 5 <= 0)
    @fact lastus(rm) --> "u $le 5"
    @addConstraint(rm, 5 <= u)
    @fact lastus(rm) --> "-u $le -5"
    @addConstraint(rm, 0 <= u - 5)
    @fact lastus(rm) --> "-u $le -5"

    @addConstraint(rm, u == sum{i*v[i], i=1:3})
    @fact lastus(rm) --> "u - v[1] - 2 v[2] - 3 v[3] $eq 0"
    @addConstraint(rm, sum{i*v[i], i=1:3} + u >= 10)
    @fact lastus(rm) --> "v[1] + 2 v[2] + 3 v[3] + u $ge 10"
end


facts("[macro] Uncertain constraints") do

    rm = RobustModel()
    @defVar(rm, x)
    @defVar(rm, y[1:5])
    @defVar(rm, z)

    @defUnc(rm, u)
    @defUnc(rm, v[1:5])
    @defUnc(rm, w)

    @addConstraint(rm, u*x <= 5)
    @fact lastuc(rm) --> "u x $le 5"

    @addConstraint(rm, x*u <= 5)
    @fact lastuc(rm) --> "u x $le 5"

    @addConstraint(rm, 5 <= u*x)
    @fact lastuc(rm) --> "-u x $le -5"

    @addConstraint(rm, 5 <= x*u)
    @fact lastuc(rm) --> "-u x $le -5"

    @addConstraint(rm, (u+w)*x + 2 + w*x <= u*z + 3)
    if VERSION > v"0.4.0-"
        @fact lastuc(rm) --> "u x + w x + w x + -u z $le 1"
    else
        @fact lastuc(rm) --> "(u + w) x + w x + -u z $le 1"
    end

    @addConstraint(rm, sum{v[i]*y[i], i=1:5; i!=3} <= 9)
    @fact lastuc(rm) --> "v[1] y[1] + v[2] y[2] + v[4] y[4] + v[5] y[5] $le 9"

    @addConstraint(rm, sum{i*(u+v[i])*(y[i]+x), i=1:2:5} <= 0)
    @fact lastuc(rm) --> "(u + v[1]) y[1] + (u + v[1]) x + (3 u + 3 v[3]) y[3] + (3 u + 3 v[3]) x + (5 u + 5 v[5]) y[5] + (5 u + 5 v[5]) x $le 0"

    foo = u*x
    @addConstraint(rm, 2 * foo <= 0)
    @fact lastuc(rm) --> "(2 u) x $le 0"
end

