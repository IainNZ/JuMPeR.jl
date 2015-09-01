#-----------------------------------------------------------------------
# JuMPeR  --  JuMP Extension for Robust Optimization
# http://github.com/IainNZ/JuMPeR.jl
#-----------------------------------------------------------------------
# Copyright (c) 2015: Iain Dunning
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
#-----------------------------------------------------------------------
# test/oracle_general_Linf.jl
# Test GeneralOracle for uncertainty sets with L∞ norms
#-----------------------------------------------------------------------

using JuMP, JuMPeR
using FactCheck
using Compat

if !(:lp_solvers in names(Main))
    println("Loading solvers...")
    include(joinpath(Pkg.dir("JuMP"),"test","solvers.jl"))
end
lp_solvers  = filter(s->(!contains(string(typeof(s)),"SCSSolver")), lp_solvers)

const TOL = 1e-4

facts("[oracle_gen_L∞] Test 1") do
for solver in lp_solvers, cuts in [true,false], flip in [true,false]
context("$(typeof(solver)), cuts=$cuts, flip=$flip") do
    m = RobustModel(solver=solver)
    @defVar(m, 0 <= x <= 10)
    @defUnc(m, 0 <= u <= 10)
    @setObjective(m, Max, 10x)
    !flip && @addConstraint(m,  u*x <=  7)
     flip && @addConstraint(m, -u*x >= -7)
    @addConstraint(m, norm(u-5, Inf) <= 2)
    solve(m, suppress_warnings=true, prefer_cuts=cuts)
    @fact getValue(x) --> roughly(1.0,TOL)
end; end; end


facts("[oracle_gen_L∞] Test 2") do
for solver in lp_solvers, cuts in [true,false],
    flip in [true,false], macr in [true,false]
context("$(typeof(solver)), cuts=$cuts, flip=$flip, macr=$macr") do
    m = RobustModel(solver=solver)
    @defVar(m, 0 <= x[i=1:5] <= 2*i)
    @defUnc(m, 0 <= u[i=1:5] <= i+4)
    @setObjective(m, Max, sum{(6-i)*x[i], i=1:5})
    !flip && @addConstraint(m,  sum{u[i]*x[i], i=1:5} <=  100)
     flip && @addConstraint(m, -sum{u[i]*x[i], i=1:5} >= -100)
    a = Float64[2, 0, 0, 2, 2];
    c = Float64[5, 0, 0, 5, 5]
    I = [1, 4, 5]
    z = convert(Vector{UAffExpr}, a.*u-c)
    !macr && @addConstraint(m, norm(z, Inf) <= 2)
     macr && @addConstraint(m, norm∞{a[i]*u[i]-c[i],i=I} <= 2)
    solve(m, suppress_warnings=true, prefer_cuts=cuts, cut_tol=1e-4)
    # max_u = [5.0, 6, 7.0, 8.0, 9.0]
    #     u = [3.5, 6, 7.0, 3.5, 3.5]
    #     x = [2.0, 4, ???, 8.0, 0.4]
    # 3.5*2.0 =  7.0 ->  7
    # 6.0*4.0 = 24.0 -> 31
    # 3.5*8.0 = 28.0 -> 59
    # x4 = 41/7
    @fact getValue(x[1]) --> roughly(2, TOL)
    @fact getValue(x[2]) --> roughly(4, TOL)
    @fact getValue(x[3]) --> roughly(41/7, TOL)
    @fact getValue(x[4]) --> roughly(8, TOL)
    @fact getValue(x[5]) --> roughly(0, TOL)
end; end; end


facts("[oracle_gen_L∞] Test 3") do
for solver in lp_solvers, cuts in [true,false],
    flip in [true,false], macr in [true,false]
context("$(typeof(solver)), cuts=$cuts, flip=$flip") do
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
    !macr && @addConstraint(m, norm(z,Inf) <= 1)
     macr && @addConstraint(m, norm∞{z[i],i=1:2} <= 1)
    solve(m, suppress_warnings=true, prefer_cuts=cuts)
    @fact getValue(x[1]) --> roughly(1.000, 1e-3)
    @fact getValue(x[2]) --> roughly(0.000, 1e-3)
end; end; end


facts("[oracle_gen_L∞] Test 4") do
for solver in lp_solvers, cuts in [true,false], flip in [true,false]
context("$(typeof(solver)), cuts=$cuts, flip=$flip") do
    m = RobustModel(solver=solver)
    @defVar(m, 0 <= x <= 10)
    @defUnc(m, 0 <= u <= 10)
    @setObjective(m, Min, 10x)
    !flip && @addConstraint(m,  x >=  u)
     flip && @addConstraint(m, -x <= -u)
    @addConstraint(m, norm(u-5,Inf) <= 2)
    solve(m, suppress_warnings=true, prefer_cuts=cuts)
    @fact getValue(x) --> roughly(7.0, TOL)
end; end; end


facts("[oracle_gen_L∞] Test 5") do
for solver in lp_solvers, cuts in [true,false], flip in [true,false]
context("$(typeof(solver)), cuts=$cuts, flip=$flip") do
    m = RobustModel(solver=solver)
    @defVar(m, 0 <= x <=  8)
    @defUnc(m, 0 <= u <= 10)

    @defVar(m, 2 <= y <= 10)
    @defUnc(m, 0 <= w <= 10)

    @setObjective(m, Max, 20x + 10y)
    !flip && @addConstraint(m,  u*x + w*y <=  10)
     flip && @addConstraint(m, -u*x - w*y >= -10)
    @addConstraint(m, norm(u - 5,Inf) <= 2)  # 3 <= u <= 7
    @addConstraint(m, norm(w - 3,Inf) <= 1)  # 2 <= w <= 4
    solve(m, suppress_warnings=true, prefer_cuts=cuts)
    @fact getValue(x) --> roughly((10-4*2)/7, TOL)
    @fact getValue(y) --> roughly(2.0, TOL)
end; end; end


facts("[oracle_gen_L∞] Test 6") do
for solver in lp_solvers, cuts in [true,false], flip in [true,false]
context("$(typeof(solver)), cuts=$cuts, flip=$flip") do
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
    @addConstraint(m, norm(u[2]-1.2,Inf) <= 0.01)
    solve(m, suppress_warnings=true, prefer_cuts=cuts)
    @fact getValue(obj) --> roughly(1.19, TOL)
end; end; end