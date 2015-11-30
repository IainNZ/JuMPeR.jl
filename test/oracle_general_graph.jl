#############################################################################
# JuMPeR
# Julia for Mathematical Programming - extension for Robust Optimization
# See http://github.com/IainNZ/JuMPeR.jl
#############################################################################
# test/oracle_general_graph.jl
# Test GeneralGraphOracle's graph functionality works
#############################################################################

using JuMP, JuMPeR
using FactCheck

if !(:lp_solvers in names(Main))
    println("Loading solvers...")
    include(joinpath(Pkg.dir("JuMP"),"test","solvers.jl"))
end
lp_solvers = filter(s->(!contains(string(typeof(s)),"SCSSolver")), lp_solvers)

facts("[oracle_gen_graph] Test 1") do
for solver in lp_solvers, cuts in [true,false]
context("$(typeof(solver)), cuts=$cuts") do
    m = RobustModel()
    r = JuMPeR.getRobust(m)
    @defUnc(m, 0 <= u[1:2,1:4] <= 1)
    for i = 1:2
        for j = 2:4
            @addConstraint(m, u[i,1] + u[i,j] >= 0)
        end
    end
    ret = JuMPeR.detect_components(r.numUncs, r.uncertaintyset)
    @fact ret[1] --> [1,1,1,1,2,2,2,2]
    @fact ret[2] --> [1,1,1,2,2,2]

    @defVar(m, x[1:2] >= 0)
    @setObjective(m, Min, sum(x))
    setDefaultOracle!(m, JuMPeR.GeneralGraphOracle())
    @addConstraint(m, x[1] >= u[1,1])
    @addConstraint(m, x[2] >= u[1,4])
    solve(m, prefer_cuts=cuts)
    @fact getValue(x[1]) --> roughly(1.0, 1e-6)
    @fact getValue(x[2]) --> roughly(1.0, 1e-6)
end; end; end


facts("[oracle_gen_graph] Test 2") do
for solver in lp_solvers, cuts in [true,false]
context("$(typeof(solver)), cuts=$cuts") do
    m = RobustModel()
    r = JuMPeR.getRobust(m)
    @defUnc(m, 2 <= u[1:3] <= 2)
    @addConstraint(m, u[1] + u[2]        >= 2)
    @addConstraint(m,        u[2] + u[3] >= 2)
    ret = JuMPeR.detect_components(r.numUncs, r.uncertaintyset)
    @fact ret[1] --> [1,1,1]
    @fact ret[2] --> [1,1]

    @defVar(m, 0 <= x[1:2] <= 10)
    @setObjective(m, Max, sum(x))
    setDefaultOracle!(m, JuMPeR.GeneralGraphOracle())
    @addConstraint(m, x[1] <= u[1])
    @addConstraint(m, x[2] <= u[3])
    solve(m, prefer_cuts=cuts)
    @fact getValue(x[1]) --> roughly(2.0, 1e-6)
    @fact getValue(x[2]) --> roughly(2.0, 1e-6)
end; end; end


facts("[oracle_gen_graph] Test 3") do
for solver in lp_solvers, cuts in [true,false]
context("$(typeof(solver)), cuts=$cuts") do
    m = RobustModel()
    r = JuMPeR.getRobust(m)
    @defUnc(m, u[1:4] >= 0)
    @addConstraint(m, u[1]        + u[3]        == 1)
    @addConstraint(m,        u[2]               == 1)
    @addConstraint(m,        u[2] + u[3]        == 1)
    @addConstraint(m,                      u[4] <= 5)
    ret = JuMPeR.detect_components(r.numUncs, r.uncertaintyset)
    @fact ret[1] --> [1,1,1,2]
    @fact ret[2] --> [1,1,1,2]

    @defVar(m, x[1:2] >= 0)
    @setObjective(m, Min, sum(x))
    setDefaultOracle!(m, JuMPeR.GeneralGraphOracle())
    @addConstraint(m, x[1] >= u[1])
    @addConstraint(m, x[2] >= u[4])
    solve(m, prefer_cuts=cuts)
    @fact getValue(x[1]) --> roughly(1.0, 1e-6)
    @fact getValue(x[2]) --> roughly(5.0, 1e-6)
end; end; end
