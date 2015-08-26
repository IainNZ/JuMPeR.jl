#############################################################################
# JuMPeR
# Julia for Mathematical Programming - extension for Robust Optimization
# See http://github.com/IainNZ/JuMPeR.jl
#############################################################################
# test/oracle_bertsim.jl
# Test BertSimOracle
#############################################################################

using JuMP, JuMPeR
using FactCheck
using Compat

if !(:lp_solvers in names(Main))
    println("Loading solvers...")
    include(joinpath(Pkg.dir("JuMP"),"test","solvers.jl"))
end
lp_solvers = filter(s->(!contains(string(typeof(s)),"SCSSolver")), lp_solvers)

facts("[oracle_bertsim] +x, +coeff") do
for solver in lp_solvers, cuts in [true]
context("$(typeof(solver)), cuts=$cuts") do
    n           = 2
    weight_low  = [1.0, 5.0]
    weight_high = [3.0, 7.0]
    values      = [0.1, 9.9]

    m = RobustModel(solver=solver)
    @fact BertSimOracle() --> not(nothing)
    setDefaultOracle!(m, BertSimOracle(1))
    
    @defVar(m, 0 <= x[1:n] <= 10)
    @defUnc(m, weight_low[i] <= u[i=1:n] <= weight_high[i])

    @setObjective(m, Max, sum{values[i] * x[i], i=1:n})

    @addConstraint(m, sum{u[i]*x[i], i=1:n} <= 21)
    
    @fact solve(m, prefer_cuts=cuts) --> :Optimal
    @fact getValue(x[1]) --> roughly(0.0,1e-6)
    @fact getValue(x[2]) --> roughly(3.0,1e-6)
end; end; end


facts("[oracle_bertsim] +x, -coeff") do
for solver in lp_solvers, cuts in [true]
context("$(typeof(solver)), cuts=$cuts") do
    n           = 2
    weight_high = [-1.0, -5.0]
    weight_low  = [-3.0, -7.0]
    values      = [0.1, 9.9]

    m = RobustModel(solver=solver)
    setDefaultOracle!(m, BertSimOracle(1))
    
    @defVar(m, 0 <= x[1:n] <= 10)
    @defUnc(m, weight_low[i] <= u[i=1:n] <= weight_high[i])

    @setObjective(m, Max, sum{ values[i] * x[i], i=1:n})

    @addConstraint(m, sum{u[i]*x[i], i=1:n} >= -21)
    
    @fact solve(m, prefer_cuts=cuts) --> :Optimal
    @fact getValue(x[1]) --> roughly(0.0,1e-6)
    @fact getValue(x[2]) --> roughly(3.0,1e-6)
end; end; end

facts("[oracle_bertsim] -x, +coeff") do
for solver in lp_solvers, cuts in [true]
context("$(typeof(solver)), cuts=$cuts") do
    n           = 2
    weight_low  = [1.0, 5.0]
    weight_high = [3.0, 7.0]
    values      = [-0.1, -9.9]

    m = RobustModel(solver=solver)
    setDefaultOracle!(m, BertSimOracle(1))
    
    @defVar(m, -10 <= x[1:n] <= 0)
    @defUnc(m, weight_low[i] <= u[i=1:n] <= weight_high[i])

    @setObjective(m, Max, sum{ values[i] * x[i], i=1:n})

    @addConstraint(m, sum{u[i]*x[i], i=1:n} >= -21)
    
    @fact solve(m, prefer_cuts=cuts) --> :Optimal
    @fact getValue(x[1]) --> roughly( 0.0,1e-6)
    @fact getValue(x[2]) --> roughly(-3.0,1e-6)
end; end; end


facts("[oracle_bertsim] -x, -coeff") do
for solver in lp_solvers, cuts in [true]
context("$(typeof(solver)), cuts=$cuts") do
    n           = 2
    weight_high = [-1.0, -5.0]
    weight_low  = [-3.0, -7.0]
    values      = [-0.1, -9.9]

    m = RobustModel(solver=solver)
    setDefaultOracle!(m, BertSimOracle(1))
    
    @defVar(m, -10 <= x[1:n] <= 0)
    @defUnc(m, weight_low[i] <= u[i=1:n] <= weight_high[i])

    @setObjective(m, Max, sum{ values[i] * x[i], i=1:n})

    @addConstraint(m, sum{u[i]*x[i], i=1:n} <= 21)
    
    @fact solve(m, prefer_cuts=cuts) --> :Optimal
    @fact getValue(x[1]) --> roughly( 0.0,1e-6)
    @fact getValue(x[2]) --> roughly(-3.0,1e-6)
end; end; end