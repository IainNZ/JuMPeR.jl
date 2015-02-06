#############################################################################
# JuMPeR
# Julia for Mathematical Programming - extension for Robust Optimization
# See http://github.com/IainNZ/JuMPeR.jl
#############################################################################
# test/operators.jl
# Testing for all operator overloads
#############################################################################

using JuMP, JuMPeR
using FactCheck

facts("[ellipse] Check the construction of the constraint data") do
    rm = RobustModel()
    @defUnc(rm, u[1:5])

    vec = [1.0*u[1]+2.0*u[2]+1.0, 5.0*u[5]+1.0*u[1]+5.0, 4.0*u[4]+4.0]
    con = build_ellipse_constraint(vec, 2.0)
    @fact con.F => [1.0 2.0 0.0 0.0;
                    1.0 0.0 5.0 0.0;
                    0.0 0.0 0.0 4.0]
    @fact con.u => [  1,  2,  5,  4]
    @fact con.g => [1.0, 5.0, 4.0]
    @fact con.Gamma => 2.0

    vec = [u[5], u[3], u[1] + 9.0]
    con = build_ellipse_constraint(vec, 1.0)
    @fact con.F => [1.0 0.0 0.0;
                    0.0 1.0 0.0;
                    0.0 0.0 1.0]
    @fact con.u => [  5,  3,  1]
    @fact con.g => [ 0., 0., 9.]
    @fact con.Gamma => 1.0

    vec = [1.0*u[1]-5]
    con = build_ellipse_constraint(vec, 3.0)
    @fact con.F => ones(1,1)
    @fact con.u =>[  1]
    @fact con.g =>[-5.]
    @fact con.Gamma => 3.0

    vec = [2.0, 0.0]
    @fact_throws build_ellipse_constraint(vec, 0.0)

    addEllipseConstraint(rm, [u[5], u[3], u[1] + 9.0], 1.0)
end

# Load some solvers
grb = isdir(Pkg.dir("Gurobi"))
eco = isdir(Pkg.dir("ECOS"))
if grb; import Gurobi; end
if eco; import ECOS; end
ell_solvers = Any[]
grb && push!(ell_solvers, Gurobi.GurobiSolver(OutputFlag=0))
eco && push!(ell_solvers, ECOS.ECOSSolver(verbose=0))

facts("[ellipse] Test using ellipsoidal sets") do
for ell_solver in ell_solvers, cuts in [true,false], flip in [true,false]
context("Testing with $(ell_solver), cuts=$cuts, flip=$flip") do
context("Test 1 - 1 var, 1 lhs unc") do
    m = RobustModel(solver=ell_solver)
    @defVar(m, 0 <= x <= 10)
    @defUnc(m, 0 <= u <= 10)
    @setObjective(m, Max, 10x)
    !flip && addConstraint(m,  u*x <=  7)
     flip && addConstraint(m, -u*x >= -7)
    addEllipseConstraint(m, [1.0*u - 5], 2)  # 5 <= u <= 7
    solveRobust(m, suppress_warnings=true, prefer_cuts=cuts)
    @fact getValue(x) => roughly(1.0,1e-6)
end
context("Test 2 - 5 var, 5 unc") do
    m = RobustModel(solver=ell_solver)
    @defVar(m, 0 <= x[i=1:5] <= 2*i)
    @defUnc(m, 0 <= u[i=1:5] <= i+4)
    @setObjective(m, Max, sum{(6-i)*x[i], i=1:5})
    !flip && addConstraint(m,  sum([ u[i]*x[i] for i=1:5 ]) <=  100)
     flip && addConstraint(m, -sum([ u[i]*x[i] for i=1:5 ]) >= -100)
    addEllipseConstraint(m, [3.0*u[1]-5, 1.0*u[5]-5, 2.0*u[4]-5], 1)
    solveRobust(m, suppress_warnings=true, prefer_cuts=cuts)
    @fact getValue(x[1]) => roughly(2.0,1e-6)
    @fact getValue(x[2]) => roughly(4.0,1e-6)
    @fact getValue(x[3]) => roughly(6.0,1e-6)
    @fact getValue(x[4]) => roughly(8.0,1e-6)
    @fact getValue(x[5]) => roughly(1.283,1e-2)
end
context("Test 3 - Polyhedral + ellipse") do
    m = RobustModel(solver=ell_solver)
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
    solveRobust(m, suppress_warnings=true, prefer_cuts=cuts)
    @fact getValue(x[1]) => roughly(1.000, 1e-3)
    @fact getValue(x[2]) => roughly(0.000, 1e-3)
end
context("Test 4 - 1 var, 1 rhs unc") do
    m = RobustModel(solver=ell_solver)
    @defVar(m, 0 <= x <= 10)
    @defUnc(m, 0 <= u <= 10)
    @setObjective(m, Min, 10x)
    !flip && addConstraint(m,  x >=  u)
     flip && addConstraint(m, -x <= -u)
    addEllipseConstraint(m, [1.0*u - 5], 2)  # 5 <= u <= 7
    solveRobust(m, suppress_warnings=true, prefer_cuts=cuts)
    @fact getValue(x) => roughly(7.0, 1e-6)
end
context("Test 5 - 2 var, 2 unc, 2 ellipse") do
    m = RobustModel(solver=ell_solver)
    @defVar(m, 0 <= x <=  8)
    @defUnc(m, 0 <= u <= 10)

    @defVar(m, 2 <= y <= 10)
    @defUnc(m, 0 <= w <= 10)

    @setObjective(m, Max, 20x + 10y)
    !flip && addConstraint(m,  u*x + w*y <=  10)
     flip && addConstraint(m, -u*x - w*y >= -10)
    addEllipseConstraint(m, [u - 5], 2)  # 5 <= u <= 7
    addEllipseConstraint(m, [w - 3], 1)  # 2 <= w <= 4
    solveRobust(m, suppress_warnings=true, prefer_cuts=cuts) #, debug_printreform=!cuts, debug_printfinal=!cuts)
    @fact getValue(x) => roughly((10-4*2)/7, 1e-5)
    @fact getValue(y) => roughly(2.0, 1e-5)
end
context("Test 6 - dfagnan provided") do
    m = RobustModel(solver=ell_solver)
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
    solveRobust(m, suppress_warnings=true, prefer_cuts=cuts)
    @fact getValue(obj) => roughly(1.19, 1e-6)
end
end # outer context
end # for
end # facts