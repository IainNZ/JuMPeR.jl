#-----------------------------------------------------------------------
# JuMPeR  --  JuMP Extension for Robust Optimization
# http://github.com/IainNZ/JuMPeR.jl
#-----------------------------------------------------------------------
# Copyright (c) 2016: Iain Dunning
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
#-----------------------------------------------------------------------
# test/uncsets_basic_Linf.jl
# Test BasicUncertaintySet for uncertainty sets with L∞ norms.
#-----------------------------------------------------------------------

using JuMP, JuMPeR
using BaseTestNext

const TOL = 1e-4

if !(:lp_solvers in names(Main))
    print_with_color(:magenta, "Loading solvers...\n")
    include(joinpath(Pkg.dir("JuMP"),"test","solvers.jl"))
end
lp_solvers  = filter(s->(!contains(string(typeof(s)),"SCSSolver")), lp_solvers)
solver_name(solver) = split(string(typeof(solver)),".")[2]

@testset "BasicUncertaintySet L∞ norm" begin
print_with_color(:yellow, "BasicUncertaintySet L∞ norm...\n")
@testset "LPs with $(solver_name(solver)), cuts=$cuts, flip=$flip" for
                        solver in lp_solvers, cuts in [true,false],
                        flip in [true,false]

    @testset "Test 1" begin
        m = RobustModel(solver=solver)
        @variable(m, 0 <= x <= 10)
        @uncertain(m, 0 <= u <= 10)
        @objective(m, Max, 10x)
        !flip && @constraint(m,  u*x <=  7)
         flip && @constraint(m, -u*x >= -7)
        @constraint(m, norm(u-5, Inf) <= 2)
        @test solve(m, suppress_warnings=true, prefer_cuts=cuts) == :Optimal
        @test isapprox(getvalue(x), 1.0, atol=TOL)
    end  # "Test 1"


    @testset "Test 2, macr=$macr" for macr in [true,false]
        m = RobustModel(solver=solver)
        @variable(m, 0 <= x[i=1:5] <= 2*i)
        @uncertain(m, 0 <= u[i=1:5] <= i+4)
        @objective(m, Max, sum{(6-i)*x[i], i=1:5})
        !flip && @constraint(m,  sum{u[i]*x[i], i=1:5} <=  100)
         flip && @constraint(m, -sum{u[i]*x[i], i=1:5} >= -100)
        a = Float64[2, 0, 0, 2, 2];
        c = Float64[5, 0, 0, 5, 5]
        I = [1, 4, 5]
        z = convert(Vector{JuMPeR.UncExpr}, a.*u-c)
        !macr && @constraint(m, norm(z, Inf) <= 2)
         macr && @constraint(m, norm∞{a[i]*u[i]-c[i],i=I} <= 2)
        solve(m, suppress_warnings=true, prefer_cuts=cuts, cut_tol=1e-4)
        # max_u = [5.0, 6, 7.0, 8.0, 9.0]
        #     u = [3.5, 6, 7.0, 3.5, 3.5]
        #     x = [2.0, 4, ???, 8.0, 0.4]
        # 3.5*2.0 =  7.0 ->  7
        # 6.0*4.0 = 24.0 -> 31
        # 3.5*8.0 = 28.0 -> 59
        # x4 = 41/7
        @test isapprox(getvalue(x[1]), 2, atol=TOL)
        @test isapprox(getvalue(x[2]), 4, atol=TOL)
        @test isapprox(getvalue(x[3]), 41/7, atol=TOL)
        @test isapprox(getvalue(x[4]), 8, atol=TOL)
        @test isapprox(getvalue(x[5]), 0, atol=TOL)
    end  # "Test 2"

    @testset "Test 3, macr=$macr" for macr in [true,false]
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
        !macr && @constraint(m, norm(z,Inf) <= 1)
         macr && @constraint(m, norm∞{z[i],i=1:2} <= 1)
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
        @constraint(m, norm(u-5,Inf) <= 2)
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
        @constraint(m, norm(u - 5,Inf) <= 2)  # 3 <= u <= 7
        @constraint(m, norm(w - 3,Inf) <= 1)  # 2 <= w <= 4
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
        @constraint(m, norm(u[2]-1.2,Inf) <= 0.01)
        @test solve(m, suppress_warnings=true, prefer_cuts=cuts) == :Optimal
        @test isapprox(getvalue(obj), 1.19, atol=TOL)
    end  # "Test 6"
end  # "LPs with"
end  # "BasicUncertaintySet L∞ norm"
