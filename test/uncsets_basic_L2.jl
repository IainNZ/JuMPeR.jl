#-----------------------------------------------------------------------
# JuMPeR  --  JuMP Extension for Robust Optimization
# http://github.com/IainNZ/JuMPeR.jl
#-----------------------------------------------------------------------
# Copyright (c) 2016: Iain Dunning
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
#-----------------------------------------------------------------------
# test/uncsets_basic_L2.jl
# Test BasicUncertaintySet for uncertainty sets with L2 norms.
#-----------------------------------------------------------------------

using JuMP, JuMPeR
using BaseTestNext

const TOL = 1e-4

if !(:lp_solvers in names(Main))
    print_with_color(:magenta, "Loading solvers...\n")
    include(joinpath(Pkg.dir("JuMP"),"test","solvers.jl"))
end
soc_solvers = filter(s->(!contains(string(typeof(s)),"SCSSolver")), soc_solvers)
solver_name(solver) = split(string(typeof(solver)),".")[2]

@testset "BasicUncertaintySet L2 norm" begin
print_with_color(:yellow, "BasicUncertaintySet L2 norm...\n")
@testset "SOCPs with $(solver_name(solver)), cuts=$cuts, flip=$flip" for
                        solver in soc_solvers, cuts in [true,false],
                        flip in [true,false]

    @testset "Test 1" begin
        m = RobustModel(solver=solver)
        @defVar(m, 0 <= x <= 10)
        @defUnc(m, 0 <= u <= 10)
        @setObjective(m, Max, 10x)
        !flip && @addConstraint(m,  u*x <=  7)
         flip && @addConstraint(m, -u*x >= -7)
        @addConstraint(m, norm(u-5) <= 2)
        @test solve(m, suppress_warnings=true, prefer_cuts=cuts) == :Optimal
        @test isapprox(getValue(x), 1.0, atol=TOL)
    end  # "Test 1"

    @testset "Test 2" begin
        m = RobustModel(solver=solver)
        @defVar(m, 0 <= x[i=1:5] <= 2*i)
        @defUnc(m, 0 <= u[i=1:5] <= i+4)
        @setObjective(m, Max, sum{(6-i)*x[i], i=1:5})
        !flip && @addConstraint(m,  sum{u[i]*x[i], i=1:5} <=  100)
         flip && @addConstraint(m, -sum{u[i]*x[i], i=1:5} >= -100)
        a = [3, 0, 0, 2, 1];
        I = [1, 5, 4]
        @addConstraint(m, norm2{a[i]*u[i]-5,i=I} <= 1)
        @test solve(m, suppress_warnings=true, prefer_cuts=cuts) == :Optimal
        @test isapprox(getValue(x[1]), 2.0, atol=TOL)
        @test isapprox(getValue(x[2]), 4.0, atol=TOL)
        @test isapprox(getValue(x[3]), 6.0, atol=TOL)
        @test isapprox(getValue(x[4]), 8.0, atol=TOL)
        @test isapprox(getValue(x[5]), 1.283, atol=TOL*100)
    end  # "Test 2"

    @testset "Test 3" begin
        m = RobustModel(solver=solver)
        @defVar(m, 0 <= x[1:2] <= 10)
        @defUnc(m, u[1:2])
        @defUnc(m, z[1:2])
        @setObjective(m, Min, 1x[1] + 2x[2])
        !flip && @addConstraint(m,  u[1]*x[1] + u[2]*x[2] >=  5)
         flip && @addConstraint(m, -u[1]*x[1] - u[2]*x[2] <= -5)
        # Uncertainty set
        @addConstraint(m, u[1] == 5.0*z[1]            + 10.0)
        @addConstraint(m, u[2] == 3.0*z[1] - 2.0*z[2] +  3.0)
        @addConstraint(m, norm(z) <= 1)
        @test solve(m, suppress_warnings=true, prefer_cuts=cuts) == :Optimal
        @test isapprox(getValue(x[1]), 1.0, atol=TOL/10)
        @test isapprox(getValue(x[2]), 0.0, atol=TOL/10)
    end  # "Test 3"

    @testset "Test 4" begin
        m = RobustModel(solver=solver)
        @defVar(m, 0 <= x <= 10)
        @defUnc(m, 0 <= u <= 10)
        @setObjective(m, Min, 10x)
        !flip && @addConstraint(m,  x >=  u)
         flip && @addConstraint(m, -x <= -u)
        @addConstraint(m, norm(u-5) <= 2)
        @test solve(m, suppress_warnings=true, prefer_cuts=cuts) == :Optimal
        @test isapprox(getValue(x), 7.0, atol=TOL)
    end  # "Test 4"

    @testset "Test 5" begin
        m = RobustModel(solver=solver)
        @defVar(m, 0 <= x <=  8)
        @defUnc(m, 0 <= u <= 10)
        @defVar(m, 2 <= y <= 10)
        @defUnc(m, 0 <= w <= 10)
        @setObjective(m, Max, 20x + 10y)
        !flip && @addConstraint(m,  u*x + w*y <=  10)
         flip && @addConstraint(m, -u*x - w*y >= -10)
        @addConstraint(m, norm(u - 5) <= 2)  # 5 <= u <= 7
        @addConstraint(m, norm(w - 3) <= 1)  # 2 <= w <= 4
        @test solve(m, suppress_warnings=true, prefer_cuts=cuts) == :Optimal
        @test isapprox(getValue(x), (10-4*2)/7, atol=TOL)
        @test isapprox(getValue(y), 2.0, atol=TOL)
    end  # "Test 5"

    @testset "Test 6" begin
        m = RobustModel(solver=solver)
        @defVar(m, 0 <= y <= 100)
        @defVar(m, 0 <= z <= 100)
        @defVar(m, -100 <= obj <= 100)
        @defUnc(m, 0 <= u[1:2] <= 100)
        @setObjective(m, Max, obj)
        if !flip
            @addConstraint(m, obj <= (z+u[2]*y))
            @addConstraint(m, z <= (1-u[1]*y))
        else
            @addConstraint(m, -obj >= -(z+u[2]*y))
            @addConstraint(m, -z >= -(1-u[1]*y))
        end
        @addConstraint(m, u[1] == 1)
        @addConstraint(m, norm(u[2]-1.2) <= 0.01)
        @test solve(m, suppress_warnings=true, prefer_cuts=cuts) == :Optimal
        @test isapprox(getValue(obj), 1.19, atol=TOL)
    end  # "Test 6"
end  # "SOCPs with"
end  # "BasicUncertaintySet L2 norm"
