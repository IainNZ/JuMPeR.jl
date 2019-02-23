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
using Test

TOL = 1e-4

if !(:lp_solvers in names(Main))
    printstyled("Loading solvers...\n", color = :magenta)
    include(joinpath(dirname(pathof(JuMP)),"..","test","solvers.jl"))
end
soc_solvers = filter(s->(!occursin("SCSSolver", string(typeof(s)))), soc_solvers)
solver_name(solver) = split(string(typeof(solver)),".")[2]

@testset "BasicUncertaintySet L2 norm" begin
printstyled("BasicUncertaintySet L2 norm...\n", color = :yellow)
@testset "SOCPs with $(solver_name(solver)), cuts=$cuts, flip=$flip" for
                        solver in soc_solvers, cuts in [true,false],
                        flip in [true,false]

    @testset "Test 1" begin
        m = RobustModel(solver=solver)
        @variable(m, 0 <= x <= 10)
        @uncertain(m, 0 <= u <= 10)
        @objective(m, Max, 10x)
        !flip && @constraint(m,  u*x <=  7)
         flip && @constraint(m, -u*x >= -7)
        @constraint(m, norm(u-5) <= 2)
        @test solve(m, suppress_warnings=true, prefer_cuts=cuts) == :Optimal
        @test isapprox(getvalue(x), 1.0, atol=TOL)
    end  # "Test 1"

    @testset "Test 2" begin
        m = RobustModel(solver=solver)
        @variable(m, 0 <= x[i=1:5] <= 2*i)
        @uncertain(m, 0 <= u[i=1:5] <= i+4)
        @objective(m, Max, sum((6-i)*x[i] for i=1:5))
        !flip && @constraint(m,  sum(u[i]*x[i] for i=1:5) <=  100)
         flip && @constraint(m, -sum(u[i]*x[i] for i=1:5) >= -100)
        a = [3, 0, 0, 2, 1];
        I = [1, 5, 4]
        @constraint(m, norm(a[i]*u[i]-5 for i=I) <= 1)
        @test solve(m, suppress_warnings=true, prefer_cuts=cuts) == :Optimal
        @test isapprox(getvalue(x[1]), 2.0, atol=TOL)
        @test isapprox(getvalue(x[2]), 4.0, atol=TOL)
        @test isapprox(getvalue(x[3]), 6.0, atol=TOL)
        @test isapprox(getvalue(x[4]), 8.0, atol=TOL)
        @test isapprox(getvalue(x[5]), 1.283, atol=TOL*100)
    end  # "Test 2"

    @testset "Test 3" begin
        m = RobustModel(solver=solver)
        @variable(m, 0 <= x[1:2] <= 10)
        @uncertain(m, u[1:2])
        @uncertain(m, z[1:2])
        @objective(m, Min, 1x[1] + 2x[2])
        !flip && @constraint(m,  u[1]*x[1] + u[2]*x[2] >=  5)
         flip && @constraint(m, -u[1]*x[1] - u[2]*x[2] <= -5)
        # Uncertainty set
        @constraint(m, u[1] == 5.0*z[1]            + 10.0)
        @constraint(m, u[2] == 3.0*z[1] - 2.0*z[2] +  3.0)
        @constraint(m, norm(z) <= 1)
        @test solve(m, suppress_warnings=true, prefer_cuts=cuts) == :Optimal
        @test isapprox(getvalue(x[1]), 1.0, atol=TOL/10)
        @test isapprox(getvalue(x[2]), 0.0, atol=TOL/10)
    end  # "Test 3"

    @testset "Test 4" begin
        m = RobustModel(solver=solver)
        @variable(m, 0 <= x <= 10)
        @uncertain(m, 0 <= u <= 10)
        @objective(m, Min, 10x)
        !flip && @constraint(m,  x >=  u)
         flip && @constraint(m, -x <= -u)
        @constraint(m, norm(u-5) <= 2)
        @test solve(m, suppress_warnings=true, prefer_cuts=cuts) == :Optimal
        @test isapprox(getvalue(x), 7.0, atol=TOL)
    end  # "Test 4"

    @testset "Test 5" begin
        m = RobustModel(solver=solver)
        @variable(m, 0 <= x <=  8)
        @uncertain(m, 0 <= u <= 10)
        @variable(m, 2 <= y <= 10)
        @uncertain(m, 0 <= w <= 10)
        @objective(m, Max, 20x + 10y)
        !flip && @constraint(m,  u*x + w*y <=  10)
         flip && @constraint(m, -u*x - w*y >= -10)
        @constraint(m, norm(u - 5) <= 2)  # 5 <= u <= 7
        @constraint(m, norm(w - 3) <= 1)  # 2 <= w <= 4
        @test solve(m, suppress_warnings=true, prefer_cuts=cuts) == :Optimal
        @test isapprox(getvalue(x), (10-4*2)/7, atol=TOL)
        @test isapprox(getvalue(y), 2.0, atol=TOL)
    end  # "Test 5"

    @testset "Test 6" begin
        m = RobustModel(solver=solver)
        @variable(m, 0 <= y <= 100)
        @variable(m, 0 <= z <= 100)
        @variable(m, -100 <= obj <= 100)
        @uncertain(m, 0 <= u[1:2] <= 100)
        @objective(m, Max, obj)
        if !flip
            @constraint(m, obj <= (z+u[2]*y))
            @constraint(m, z <= (1-u[1]*y))
        else
            @constraint(m, -obj >= -(z+u[2]*y))
            @constraint(m, -z >= -(1-u[1]*y))
        end
        @constraint(m, u[1] == 1)
        @constraint(m, norm(u[2]-1.2) <= 0.01)
        @test solve(m, suppress_warnings=true, prefer_cuts=cuts) == :Optimal
        @test isapprox(getvalue(obj), 1.19, atol=TOL)
    end  # "Test 6"
end  # "SOCPs with"
end  # "BasicUncertaintySet L2 norm"
