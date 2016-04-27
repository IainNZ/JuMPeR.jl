#-----------------------------------------------------------------------
# JuMPeR  --  JuMP Extension for Robust Optimization
# http://github.com/IainNZ/JuMPeR.jl
#-----------------------------------------------------------------------
# Copyright (c) 2016: Iain Dunning
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
#-----------------------------------------------------------------------
# test/macro.jl
# Testing JuMP-defined macros still work with Uncertains etc.
#-----------------------------------------------------------------------

using JuMP, JuMPeR
using BaseTestNext

lastuc(rm) = JuMP.con_str(JuMP.REPLMode, JuMPeR.get_robust(rm).unc_constraints[end])
lastus(rm) = JuMP.con_str(JuMP.REPLMode, JuMPeR.get_robust(rm).default_uncset.linear_constraints[end])
le, eq, ge = JuMP.repl[:leq], JuMP.repl[:eq], JuMP.repl[:geq]

@testset "Macros" begin
print_with_color(:yellow, "Macros...\n")

@testset "Uncertain parameters" begin
    rm = RobustModel()
    @uncertain(rm, 5 >= u >= 1)
    @test JuMPeR.get_robust(rm).unc_lower[1] == 1
    @test JuMPeR.get_robust(rm).unc_upper[1] == 5
    @test string(zero(u)) == "0"
    @test string(one(u)) == "1"
    @test sprint(print,rm) == """Min 0
Subject to
Uncertain constraints:
Uncertain parameters:
1 ≤ u ≤ 5
"""
end

@testset "Unc set constraints" begin
    rm = RobustModel()
    @uncertain(rm, u)
    @uncertain(rm, v[1:3])

    @constraint(rm, u <= 5)
    @test lastus(rm) == "u $le 5"
    @constraint(rm, u - 5 <= 0)
    @test lastus(rm) == "u $le 5"
    @constraint(rm, 5 <= u)
    @test lastus(rm) == "-u $le -5"
    @constraint(rm, 0 <= u - 5)
    @test lastus(rm) == "-u $le -5"

    @constraint(rm, u == sum{i*v[i], i=1:3})
    @test lastus(rm) == "u - v[1] - 2 v[2] - 3 v[3] $eq 0"
    @constraint(rm, sum{i*v[i], i=1:3} + u >= 10)
    @test lastus(rm) == "v[1] + 2 v[2] + 3 v[3] + u $ge 10"
end


@testset "Uncertain constraints" begin
    rm = RobustModel()
    @variable(rm, x)
    @variable(rm, y[1:5])
    @variable(rm, z)

    @uncertain(rm, u)
    @uncertain(rm, v[1:5])
    @uncertain(rm, w)

    @constraint(rm, u*x <= 5)
    @test lastuc(rm) == "u x $le 5"

    @constraint(rm, x*u <= 5)
    @test lastuc(rm) == "u x $le 5"

    @constraint(rm, 5 <= u*x)
    @test lastuc(rm) == "-u x $le -5"

    @constraint(rm, 5 <= x*u)
    @test lastuc(rm) == "-u x $le -5"

    @constraint(rm, (u+w)*x + 2 + w*x <= u*z + 3)
    @test lastuc(rm) == "u x + w x + w x + -u z $le 1"

    @constraint(rm, sum{v[i]*y[i], i=1:5; i!=3} <= 9)
    @test lastuc(rm) == "v[1] y[1] + v[2] y[2] + v[4] y[4] + v[5] y[5] $le 9"

    @constraint(rm, sum{i*(u+v[i])*(y[i]+x), i=1:2:5} <= 0)
    @test lastuc(rm) == "(u + v[1]) y[1] + (u + v[1]) x + (3 u + 3 v[3]) y[3] + (3 u + 3 v[3]) x + (5 u + 5 v[5]) y[5] + (5 u + 5 v[5]) x $le 0"

    foo = u*x
    @constraint(rm, 2 * foo <= 0)
    @test lastuc(rm) == "(2 u) x $le 0"
end

end  # "Macros"
