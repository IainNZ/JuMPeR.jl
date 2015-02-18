#############################################################################
# JuMPeR
# Julia for Mathematical Programming - extension for Robust Optimization
# See http://github.com/IainNZ/JuMPeR.jl
#############################################################################
# test/oracle_general.jl
# Test GeneralOracle for polyhedral and ellipsoidal uncertainty sets
#############################################################################

using JuMP, JuMPeR
using FactCheck

if !(:lp_solvers in names(Main))
    println("Loading solvers...")
    include(joinpath(Pkg.dir("JuMP"),"test","solvers.jl"))
end

facts("[oracle_gen_poly] Test 1") do
for solver in lp_solvers, cuts in [true,false]
context("$(typeof(solver)), cuts=$cuts") do
    m = RobustModel(solver=solver)
    @defVar(m, x[1:2] >= 0)
    @defUnc(m, 0.3 <= u <= 0.5)
    @setObjective(m, Max, x[1] + x[2])
    @addConstraint(m, u*x[1] + 1*x[2] <= 2)
    @addConstraint(m, 1*x[1] + 1*x[2] <= 6)
    @fact solve(m, prefer_cuts=cuts) => :Optimal
    @fact getValue(x[1]) => roughly(4.0,1e-6)
    @fact getValue(x[2]) => roughly(0.0,1e-6)
end; end; end


facts("[oracle_gen_poly] Test 2") do
for solver in lp_solvers, cuts in [true,false]
context("$(typeof(solver)), cuts=$cuts") do
    m = RobustModel(solver=solver)
    @defVar(m, x[1:2] >= 0)
    @defUnc(m, 0.3 <= u1 <= 0.5)
    @defUnc(m, 0.0 <= u2 <= 2.0)
    @setObjective(m, Max, x[1] + x[2])
    @addConstraint(m, u1*x[1] + 1*x[2] <= 2)
    @addConstraint(m, u2*x[1] + 1*x[2] <= 6)
    @fact solve(m, prefer_cuts=cuts, add_box=cuts?1e2:false) => :Optimal
    @fact getValue(x[1]) => roughly(2.0+2.0/3.0,1e-6)
    @fact getValue(x[2]) => roughly(    2.0/3.0,1e-6)
end; end; end


facts("[oracle_gen_poly] Test 2-IP") do
for solver in lazy_solvers, cuts in [true,false]
context("$(typeof(solver)), cuts=$cuts") do
    m = RobustModel(solver=solver)
    @defVar(m, x[1:2] >= 0, Int)
    @defUnc(m, 0.3 <= u1 <= 0.5)
    @defUnc(m, 0.0 <= u2 <= 2.0)
    @setObjective(m, Max, 1.1*x[1] + x[2])
    @addConstraint(m, u1*x[1] + 1*x[2] <= 2)
    @addConstraint(m, u2*x[1] + 1*x[2] <= 6)
    @fact solve(m, prefer_cuts=cuts, add_box=cuts?1e2:false) => :Optimal
    @fact getValue(x[1]) => roughly(3.0,1e-6)
    @fact getValue(x[2]) => roughly(0.0,1e-6)
end; end; end


facts("[oracle_gen_poly] Test 3") do
for solver in lp_solvers, cuts in [true,false]
context("$(typeof(solver)), cuts=$cuts") do
    m = RobustModel(solver=solver)
    @defVar(m, x[1:2] >= 0)
    @defUnc(m, 0.3 <= u1 <= 1.5)
    @defUnc(m, 0.5 <= u2 <= 1.5)
    @setObjective(m, Max, x[1] + x[2])
    @addConstraint(m, u1*x[1] <= 3)
    @addConstraint(m, u2*x[2] <= 1)
    @addConstraint(m, (2.0*u1-2.0) + (4.0*u2-2.0) <= +1)
    @addConstraint(m, (2.0*u1-2.0) + (4.0*u2-2.0) >= -1)
    @fact solve(m, prefer_cuts=cuts, add_box=cuts?1e2:false) => :Optimal
    @fact getValue(x[1]) => roughly(2.0,1e-6)
    @fact getValue(x[2]) => roughly(10.0/11.0,1e-6)
end; end; end


facts("[oracle_gen_poly] Test 3-IP") do
for solver in lazy_solvers, cuts in [true,false]
context("$(typeof(solver)), cuts=$cuts") do
    m = RobustModel(solver=solver)
    @defVar(m, x[1:2] >= 0, Int)
    @defUnc(m, 0.3 <= u1 <= 1.5)
    @defUnc(m, 0.5 <= u2 <= 1.5)
    @setObjective(m, Max, x[1] + x[2])
    @addConstraint(m, u1*x[1] <= 3)
    @addConstraint(m, u2*x[2] <= 1)
    @addConstraint(m, (2.0*u1-2.0) + (4.0*u2-2.0) <= +1)
    @addConstraint(m, (2.0*u1-2.0) + (4.0*u2-2.0) >= -1)
    @fact solve(m, prefer_cuts=cuts, add_box=cuts?1e2:false) => :Optimal
    @fact getValue(x[1]) => roughly(2.0,1e-6)
    @fact getValue(x[2]) => roughly(0.0,1e-6)
end; end; end


facts("[oracle_gen_poly] Test 4") do
for solver in lp_solvers, cuts in [true,false]
context("$(typeof(solver)), cuts=$cuts") do
    m = RobustModel(solver=solver)
    @defVar(m, x >= 0)
    @defUnc(m, u <= 4.0)
    @addConstraint(m, u >= 3.0)
    @setObjective(m, Max, 1.0x)
    @addConstraint(m, x <= u)
    @fact solve(m, prefer_cuts=cuts, add_box=cuts?1e2:false) => :Optimal
    @fact getValue(x) => roughly(3.0,1e-6)
end; end; end


facts("[oracle_gen_poly] Test 5") do
for solver in lp_solvers, cuts in [true,false]
context("$(typeof(solver)), cuts=$cuts") do
    m = RobustModel(solver=solver)
    @defUnc(m, u[1:5] >= 0)
    for ix = 1:5
        @addConstraint(m, u[ix] <= float(ix))
    end
    @addConstraint(m, sum(u) <= 5)
    @defVar(m, t >= 0)
    @addConstraint(m, t >= sum(u))
    @setObjective(m, Min, t)
    @fact solve(m, prefer_cuts=cuts) => :Optimal
    @fact getObjectiveValue(m) => roughly(5.0,1e-6)
end; end; end


facts("[oracle_gen_poly] Test 6") do
for solver in lp_solvers, cuts in [true,false], variant in 0:7
context("$(typeof(solver)), cuts=$cuts, variant=$variant") do
    rm = RobustModel(solver=solver)
    @defUnc(rm, u >=0)
    @addConstraint(rm, u <=0)
    @defVar(rm, x >=0)
    @defVar(rm, shed >=0)
    @setObjective(rm, Min, x + shed)
    variant == 0 && @addConstraint(rm, x - u + 3.46 - shed <= 0)
    variant == 1 && @addConstraint(rm, x - u + 3.46 <= shed)
    variant == 2 && @addConstraint(rm, x - u <= shed - 3.46)
    variant == 3 && @addConstraint(rm, x <= u + shed - 3.46)
    variant == 4 && @addConstraint(rm, 0 <= -x + u + shed - 3.46)
    variant == 5 && @addConstraint(rm, 3.46 <= -x + u + shed)
    variant == 6 && @addConstraint(rm, 3.46 - shed <= -x + u)
    variant == 7 && @addConstraint(rm, 3.46 + x <= shed + u)
    @fact solve(rm, prefer_cuts=cuts) => :Optimal
    @fact getObjectiveValue(rm) => roughly(3.46,1e-6)
end; end; end


facts("[oracle_gen_poly] Test 7 (integer uncertainty set)") do
for solver in lazy_solvers
context("$(typeof(solver))") do
    rm = RobustModel(solver=solver)
    @defVar(rm, 0 <= x <= 2)
    @defUnc(rm, 0.5 <= u <= 1.5, Int)
    @setObjective(rm, Max, x)
    @addConstraint(rm, x <= u)
    @fact solve(rm, prefer_cuts=true) => :Optimal
    @fact getValue(x) => roughly(1.0,1e-6)
end; end; end


facts("[oracle_gen_poly] Test 8") do
for solver in lp_solvers, cuts in [true,false]
context("$(typeof(solver)), cuts=$cuts") do
    m = RobustModel(solver=solver)
    @defVar(m, x >= 0)
    @defUnc(m, 0.5 <= u <= 0.5)
    @setObjective(m, Max, x)
    @addConstraint(m, u*x + u <= 2)
    @fact solve(m, prefer_cuts=cuts, add_box=cuts?1e2:false) => :Optimal
    @fact getValue(x) => roughly(3.0,1e-6)
end; end; end


facts("[oracle_gen_ell] Test 1") do
for solver in soc_solvers, cuts in [true,false], flip in [true,false]
context("$(typeof(solver)), cuts=$cuts, flip=$flip") do
    m = RobustModel(solver=solver)
    @defVar(m, 0 <= x <= 10)
    @defUnc(m, 0 <= u <= 10)
    @setObjective(m, Max, 10x)
    !flip && addConstraint(m,  u*x <=  7)
     flip && addConstraint(m, -u*x >= -7)
    addEllipseConstraint(m, [1.0*u - 5], 2)  # 5 <= u <= 7
    solve(m, suppress_warnings=true, prefer_cuts=cuts)
    @fact getValue(x) => roughly(1.0,1e-6)
end; end; end


facts("[oracle_gen_ell] Test 2") do
for solver in soc_solvers, cuts in [true,false], flip in [true,false]
context("$(typeof(solver)), cuts=$cuts, flip=$flip") do
    m = RobustModel(solver=solver)
    @defVar(m, 0 <= x[i=1:5] <= 2*i)
    @defUnc(m, 0 <= u[i=1:5] <= i+4)
    @setObjective(m, Max, sum{(6-i)*x[i], i=1:5})
    !flip && addConstraint(m,  sum([ u[i]*x[i] for i=1:5 ]) <=  100)
     flip && addConstraint(m, -sum([ u[i]*x[i] for i=1:5 ]) >= -100)
    addEllipseConstraint(m, [3.0*u[1]-5, 1.0*u[5]-5, 2.0*u[4]-5], 1)
    solve(m, suppress_warnings=true, prefer_cuts=cuts)
    @fact getValue(x[1]) => roughly(2.0,1e-6)
    @fact getValue(x[2]) => roughly(4.0,1e-6)
    @fact getValue(x[3]) => roughly(6.0,1e-6)
    @fact getValue(x[4]) => roughly(8.0,1e-6)
    @fact getValue(x[5]) => roughly(1.283,1e-2)
end; end; end


facts("[oracle_gen_ell] Test 3") do
for solver in soc_solvers, cuts in [true,false], flip in [true,false]
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
    addEllipseConstraint(m, [z[1],z[2]], 1)
    solve(m, suppress_warnings=true, prefer_cuts=cuts)
    @fact getValue(x[1]) => roughly(1.000, 1e-3)
    @fact getValue(x[2]) => roughly(0.000, 1e-3)
end; end; end;


facts("[oracle_gen_ell] Test 4") do
for solver in soc_solvers, cuts in [true,false], flip in [true,false]
context("$(typeof(solver)), cuts=$cuts, flip=$flip") do
    m = RobustModel(solver=solver)
    @defVar(m, 0 <= x <= 10)
    @defUnc(m, 0 <= u <= 10)
    @setObjective(m, Min, 10x)
    !flip && addConstraint(m,  x >=  u)
     flip && addConstraint(m, -x <= -u)
    addEllipseConstraint(m, [1.0*u - 5], 2)  # 5 <= u <= 7
    solve(m, suppress_warnings=true, prefer_cuts=cuts)
    @fact getValue(x) => roughly(7.0, 1e-6)
end; end; end


facts("[oracle_gen_ell] Test 5") do
for solver in soc_solvers, cuts in [true,false], flip in [true,false]
context("$(typeof(solver)), cuts=$cuts, flip=$flip") do
    m = RobustModel(solver=solver)
    @defVar(m, 0 <= x <=  8)
    @defUnc(m, 0 <= u <= 10)

    @defVar(m, 2 <= y <= 10)
    @defUnc(m, 0 <= w <= 10)

    @setObjective(m, Max, 20x + 10y)
    !flip && addConstraint(m,  u*x + w*y <=  10)
     flip && addConstraint(m, -u*x - w*y >= -10)
    addEllipseConstraint(m, [u - 5], 2)  # 5 <= u <= 7
    addEllipseConstraint(m, [w - 3], 1)  # 2 <= w <= 4
    solve(m, suppress_warnings=true, prefer_cuts=cuts)
    @fact getValue(x) => roughly((10-4*2)/7, 1e-5)
    @fact getValue(y) => roughly(2.0, 1e-5)
end; end; end


facts("[oracle_gen_ell] Test 6") do
for solver in soc_solvers, cuts in [true,false], flip in [true,false]
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
    addEllipseConstraint(m,[u[2]-1.2],0.01)
    solve(m, suppress_warnings=true, prefer_cuts=cuts)
    @fact getValue(obj) => roughly(1.19, 1e-6)
end; end; end


facts("[oracle_gen] Resolve") do
for solver in lp_solvers, cuts in [true,false]
context("$(typeof(solver)), cuts=$cuts") do
    m = RobustModel(solver=solver)
    @defVar(m, B[1:2] >= 0)
    @defVar(m, S[1:2] >= 0)
    
    # Uncertainty set
    @defUnc(m, D[1:2])
    @addConstraint(m,  (D[1] - 20)/10 + (D[2] - 10)/5 <= 1.5)
    @addConstraint(m,  (D[1] - 20)/10 - (D[2] - 10)/5 <= 1.5)
    @addConstraint(m, -(D[1] - 20)/10 + (D[2] - 10)/5 <= 1.5)
    @addConstraint(m, -(D[1] - 20)/10 - (D[2] - 10)/5 <= 1.5)

    # Constraints
    for i in 1:2
        @addConstraint(m, S[i] <= D[i])
        @addConstraint(m, S[i] <= B[i])
    end

    # Objective
    @setObjective(m, Max, 3*sum(S) - sum(B))

    # First solve
    @fact solve(m, prefer_cuts=cuts, add_box=cuts?1e2:false) => :Optimal
    @fact getValue(B[1]) => roughly(5.0,1e-6)
    @fact getValue(B[2]) => roughly(2.5,1e-6)

    # Solve again with no changes
    @fact solve(m, prefer_cuts=cuts, add_box=cuts?1e2:false) => :Optimal
    @fact getValue(B[1]) => roughly(5.0,1e-6)
    @fact getValue(B[2]) => roughly(2.5,1e-6)

    # Tighten uncertainty set 
    @addConstraint(m,  (D[1] - 20)/10 + (D[2] - 10)/5 <= 1.0)
    @addConstraint(m,  (D[1] - 20)/10 - (D[2] - 10)/5 <= 1.0)
    @addConstraint(m, -(D[1] - 20)/10 + (D[2] - 10)/5 <= 1.0)
    @addConstraint(m, -(D[1] - 20)/10 - (D[2] - 10)/5 <= 1.0)
    @fact solve(m, prefer_cuts=cuts, add_box=cuts?1e2:false) => :Optimal
    @fact getValue(B[1]) => roughly(10.0,1e-6)
    @fact getValue(B[2]) => roughly( 5.0,1e-6)

    # Add a certain constraint
    @addConstraint(m, B[1] <= 8)
    @fact solve(m, prefer_cuts=cuts, add_box=cuts?1e2:false) => :Optimal
    @fact getValue(B[1]) => roughly(8.0,1e-6)
    @fact getValue(B[2]) => roughly(5.0,1e-6)

    # Add an uncertain constraint (and disambiguate objective)
    @addConstraint(m, B[1] + B[2] <= (D[1] + D[2])/2)
    @setObjective(m, Max, 3.1*S[2] + 3.0*S[1] - sum(B))
    @fact solve(m, prefer_cuts=cuts, add_box=cuts?1e2:false) => :Optimal
    @fact getValue(B[1]) => roughly(5.0,1e-6)
    @fact getValue(B[2]) => roughly(5.0,1e-6)
end; end; end