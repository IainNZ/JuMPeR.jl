#############################################################################
# JuMPeR
# Julia for Mathematical Programming - extension for Robust Optimization
# See http://github.com/IainNZ/JuMPeR.jl
#############################################################################
# test/scenario.jl
# Test active scenario code (still undocumented)
#############################################################################

using JuMP, JuMPeR
using FactCheck
using Compat

if !(:lp_solvers in names(Main))
    println("Loading solvers...")
    include(joinpath(Pkg.dir("JuMP"),"test","solvers.jl"))
end

facts("[scenario] Test providing scenarios") do
for solver in lp_solvers, cuts in [true,false]
context("$(typeof(solver)), cuts=$cuts") do
    rm = RobustModel(solver=solver)
    @defVar(rm, x >= 0)
    @defUnc(rm, 3 <= u <= 5)
    @defUnc(rm, 1 <= v <= 1)
    @setObjective(rm, Max, x)
    @addConstraint(rm, v*x <= u)

    # Now we get tricky and add a scenario that is actually
    # outside the uncertainty set - but we don't check that
    addScenario(rm, @compat Dict(u => 1.5, v => 3.0))
    
    # This scenario shouldn't be added because it isn't
    # able to complete define a constraint
    addScenario(rm, @compat Dict(v => 9999.0))

    @fact solve(rm, prefer_cuts=cuts) --> :Optimal
    @fact getValue(x) --> roughly(0.5, 1e-6)
end; end; end

facts("[scenario] Test retrieval of active cuts") do
for solver in lp_solvers
context("$(typeof(solver))") do
    rm = RobustModel(solver=solver)
    @defVar(rm, 0 <= x <= 100)
    @defUnc(rm, 3 <= u <= 5)
    @defUnc(rm, 1 <= v <= 1)
    @setObjective(rm, Max, x)
    mycon      = @addConstraint(rm, v*x <= u)
    myloosecon = @addConstraint(rm, u*x <= 10000)
    solve(rm, prefer_cuts=true, active_cuts=true)
    ascen = getScenario(mycon)
    @fact getUncValue(ascen, u) --> roughly(3.0, 1e-6)
    @fact getUncValue(ascen, v) --> roughly(1.0, 1e-6)
    @fact isBinding(ascen) --> true
    bscen = getScenario(myloosecon)
    @fact getUncValue(bscen, u) --> roughly(5.0, 1e-6)
    @fact getUncValue(bscen, v) --> roughly(1.0, 1e-6)
    @fact isBinding(bscen) --> false
end; end; end