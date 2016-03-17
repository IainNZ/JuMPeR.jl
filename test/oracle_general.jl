#-----------------------------------------------------------------------
# JuMPeR  --  JuMP Extension for Robust Optimization
# http://github.com/IainNZ/JuMPeR.jl
#-----------------------------------------------------------------------
# Copyright (c) 2016: Iain Dunning
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
#-----------------------------------------------------------------------
# test/oracle_general.jl
# Test GeneralOracle for general polyhedral uncertainty sets
#-----------------------------------------------------------------------

using JuMP, JuMPeR
using BaseTestNext

const TOL = 1e-4

if !(:lp_solvers in names(Main))
    print_with_color(:magenta, "Loading solvers...\n")
    include(joinpath(Pkg.dir("JuMP"),"test","solvers.jl"))
end
lp_solvers  = filter(s->(!contains(string(typeof(s)),"SCSSolver")), lp_solvers)
soc_solvers = filter(s->(!contains(string(typeof(s)),"SCSSolver")), soc_solvers)
solver_name(solver) = split(string(typeof(solver)),".")[2]

@testset "GeneralOracle polyhedral" begin
print_with_color(:yellow, "GeneralOracle polyhedral...\n")
print_with_color(:yellow, "  LP tests...\n")
@testset "LPs with $(solver_name(solver)), cuts=$cuts" for
                        solver in lp_solvers, cuts in [true,false]

    @testset "Test 1" begin
        m = RobustModel(solver=solver)
        @defVar(m, x[1:2] >= 0)
        @defUnc(m, 0.3 <= u <= 0.5)
        @setObjective(m, Max, x[1] + x[2])
        @addConstraint(m, u*x[1] + 1*x[2] <= 2)
        @addConstraint(m, 1*x[1] + 1*x[2] <= 6)
        @test solve(m, prefer_cuts=cuts) == :Optimal
        @test isapprox(getValue(x[1]), 4.0, atol=TOL)
        @test isapprox(getValue(x[2]), 0.0, atol=TOL)
    end  # "Test 1"

    @testset "Test 2" begin
        m = RobustModel(solver=solver)
        @defVar(m, x[1:2] >= 0)
        @defUnc(m, 0.3 <= u1 <= 0.5)
        @defUnc(m, 0.0 <= u2 <= 2.0)
        @setObjective(m, Max, x[1] + x[2])
        @addConstraint(m, u1*x[1] + 1*x[2] <= 2)
        @addConstraint(m, u2*x[1] + 1*x[2] <= 6)
        @test solve(m, prefer_cuts=cuts, add_box=cuts?1e2:false) == :Optimal
        @test isapprox(getValue(x[1]), 2.0+2.0/3.0, atol=TOL)
        @test isapprox(getValue(x[2]),     2.0/3.0, atol=TOL)
    end  # "Test 2"

    @testset "Test 3" begin
        m = RobustModel(solver=solver)
        @defVar(m, x[1:2] >= 0)
        @defUnc(m, 0.3 <= u1 <= 1.5)
        @defUnc(m, 0.5 <= u2 <= 1.5)
        @setObjective(m, Max, x[1] + x[2])
        @addConstraint(m, u1*x[1] <= 3)
        @addConstraint(m, u2*x[2] <= 1)
        @addConstraint(m, (2.0*u1-2.0) + (4.0*u2-2.0) <= +1)
        @addConstraint(m, (2.0*u1-2.0) + (4.0*u2-2.0) >= -1)
        @test solve(m, prefer_cuts=cuts, add_box=cuts?1e2:false) == :Optimal
        @test isapprox(getValue(x[1]), 2.0,       atol=TOL)
        @test isapprox(getValue(x[2]), 10.0/11.0, atol=TOL)
    end  # "Test 3"

    @testset "Test 4" begin
        m = RobustModel(solver=solver)
        @defVar(m, x >= 0)
        @defUnc(m, u <= 4.0)
        @addConstraint(m, u >= 3.0)
        @setObjective(m, Max, 1.0x)
        @addConstraint(m, x <= u)
        @test solve(m, prefer_cuts=cuts, add_box=cuts?1e2:false) == :Optimal
        @test isapprox(getValue(x), 3.0, atol=TOL)
    end  # "Test 4"

    @testset "Test 5" begin
        m = RobustModel(solver=solver)
        @defUnc(m, u[1:5] >= 0)
        for ix = 1:5
            @addConstraint(m, u[ix] <= float(ix))
        end
        @addConstraint(m, sum(u) <= 5)
        @defVar(m, t >= 0)
        @addConstraint(m, t >= sum(u))
        @setObjective(m, Min, t)
        @test solve(m, prefer_cuts=cuts) == :Optimal
        @test isapprox(getObjectiveValue(m), 5.0, atol=TOL)
    end  # "Test 5"

    @testset "Test 6 variant $variant" for variant in 0:7
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
        @test solve(rm, prefer_cuts=cuts) == :Optimal
        @test isapprox(getObjectiveValue(rm), 3.46, atol=TOL)
    end  # "Test 6"

    @testset "Test 7" begin
        m = RobustModel(solver=solver)
        @defVar(m, x >= 0)
        @defUnc(m, 0.5 <= u <= 0.5)
        @setObjective(m, Max, x)
        @addConstraint(m, u*x + u <= 2)
        @test solve(m, prefer_cuts=cuts, add_box=cuts?1e2:false) == :Optimal
        @test isapprox(getValue(x), 3.0, atol=TOL)
    end  # "Test 7"

    @testset "Test 8 (unbounded LP)" begin
        m = RobustModel(solver=solver)
        @defVar(m, x >= 0)
        @defUnc(m, 8 <= u <= 9)
        @setObjective(m, Max, x)
        @addConstraint(m, u*x >= 25)
        @test solve(m, prefer_cuts=cuts, suppress_warnings=true) == :Unbounded
    end  # "Test 8"

    if !contains("$(typeof(solver))","IpoptSolver")  # reports UserLimit
    @testset "Test 9 (infeasible LP)" begin
        m = RobustModel(solver=solver)
        @defVar(m, x)
        @defUnc(m, 8 <= u <= 9)
        @setObjective(m, Max, x)
        @addConstraint(m, x == 3)
        @addConstraint(m, u*x <= 25)
        @test solve(m, prefer_cuts=cuts, suppress_warnings=true) == :Infeasible
    end; end  # "Test 9"

    @testset "Test 10 (empty unc. set)" begin
        m = RobustModel(solver=solver)
        @defVar(m, x >= 0)
        @defUnc(m, 8 <= u <= 7)
        @setObjective(m, Min, x)
        @addConstraint(m, u*x >= 25)
        if cuts
            # Cutting plane problem is infeasible, so an error will be
            # thrown that we catch here.
            @test_throws Exception solve(m, prefer_cuts=cuts)
        else
            # Reformulation doesn't care, although it is a bit weird
            @test solve(m, prefer_cuts=cuts) == :Optimal
        end
    end

    if !contains("$(typeof(solver))","IpoptSolver")  # reports UserLimit
    @testset "Test 11 (unbounded unc. set)" begin
        m = RobustModel(solver=solver)
        @defVar(m, x >= 1)
        @defUnc(m, u <= 7)
        @setObjective(m, Min, x)
        @addConstraint(m, u*x >= 25)
        if cuts
            # Cutting plane problem is unbounded, so an error will be
            # thrown that we catch here
            @test_throws Exception solve(m, prefer_cuts=cuts)
        else
            @test solve(m, prefer_cuts=false, suppress_warnings=true) == :Infeasible
        end
    end; end

end  # "LPs with ..."


print_with_color(:yellow, "  MILP tests...\n")
@testset "MILPs with $(solver_name(solver)), cuts=$cuts" for
                        solver in lazy_solvers, cuts in [true,false]

    @testset "Test 1" begin
        m = RobustModel(solver=solver)
        @defVar(m, x[1:2] >= 0, Int)
        @defUnc(m, 0.3 <= u1 <= 0.5)
        @defUnc(m, 0.0 <= u2 <= 2.0)
        @setObjective(m, Max, 1.1*x[1] + x[2])
        @addConstraint(m, u1*x[1] + 1*x[2] <= 2)
        @addConstraint(m, u2*x[1] + 1*x[2] <= 6)
        @test solve(m, prefer_cuts=cuts, add_box=cuts?1e2:false) == :Optimal
        @test isapprox(getValue(x[1]), 3.0, atol=TOL)
        @test isapprox(getValue(x[2]), 0.0, atol=TOL)
    end

    @testset "Test 2" begin
        m = RobustModel(solver=solver)
        @defVar(m, x[1:2] >= 0, Int)
        @defUnc(m, 0.3 <= u1 <= 1.5)
        @defUnc(m, 0.5 <= u2 <= 1.5)
        @setObjective(m, Max, x[1] + x[2])
        @addConstraint(m, u1*x[1] <= 3)
        @addConstraint(m, u2*x[2] <= 1)
        @addConstraint(m, (2*u1-2) + (4*u2-2) <= +1)
        @addConstraint(m, (2*u1-2) + (4*u2-2) >= -1)
        @test solve(m, prefer_cuts=cuts, add_box=cuts?1e2:false) == :Optimal
        @test isapprox(getValue(x[1]), 2.0, atol=TOL)
        @test isapprox(getValue(x[2]), 0.0, atol=TOL)
    end

    if cuts
    @testset "Test 3 (integer unc. set)" begin
        rm = RobustModel(solver=solver)
        @defVar(rm, 0 <= x <= 2)
        @defUnc(rm, 0.5 <= u <= 1.5, Int)
        @setObjective(rm, Max, x)
        @addConstraint(rm, x <= u)
        @test solve(rm, prefer_cuts=true) == :Optimal
        @test isapprox(getValue(x), 1.0, atol=TOL)
    end; end

    @testset "Test 4 (infeasible MILP)" begin
        m = RobustModel(solver=solver)
        @defVar(m, x, Int)
        @defUnc(m, 8 <= u <= 9)
        @setObjective(m, Max, x)
        @addConstraint(m, x == 3)
        @addConstraint(m, u*x <= 25)  # u = 9, 9*3 = 27 > 25
        @test solve(m, prefer_cuts=cuts, suppress_warnings=true) == :Infeasible
    end

    if !contains("$(typeof(solver))","GLPK")  # reports Error
    @testset "Test 5 (unbounded MILP)" begin
        m = RobustModel(solver=solver)
        @defVar(m, x >= 0, Int)
        @defUnc(m, 8 <= u <= 9)
        @setObjective(m, Max, x)
        @addConstraint(m, u*x >= 25)  # u = 8, x >= 0, > 25
        @test solve(m, prefer_cuts=cuts, suppress_warnings=true) == :Unbounded
    end; end

end  # "MILPs with..."


@testset "Resolving with $(solver_name(solver)), cuts=$cuts" for
                        solver in lazy_solvers, cuts in [true,false]
    # solve() with RobustModels is intended to be an almost-pure operation
    # The goal of these tests is to make sure nothing weird happens.
    m = RobustModel(solver=solver)
    @defVar(m, B[1:2] >= 0)
    @defVar(m, S[1:2] >= 0)

    # Uncertainty set (manually expanded budget-type set)
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

    # Objective function
    @setObjective(m, Max, 3*sum(S) - sum(B))

    # First solve
    @test solve(m, prefer_cuts=cuts, add_box=cuts?1e2:false) == :Optimal
    @test isapprox(getValue(B[1]), 5.0, atol=TOL)
    @test isapprox(getValue(B[2]), 2.5, atol=TOL)

    # Solve again with no changes to the model
    @test solve(m, prefer_cuts=cuts, add_box=cuts?1e2:false) == :Optimal
    @test isapprox(getValue(B[1]), 5.0, atol=TOL)
    @test isapprox(getValue(B[2]), 2.5, atol=TOL)

    # Tighten uncertainty set
    @addConstraint(m,  (D[1] - 20)/10 + (D[2] - 10)/5 <= 1)
    @addConstraint(m,  (D[1] - 20)/10 - (D[2] - 10)/5 <= 1)
    @addConstraint(m, -(D[1] - 20)/10 + (D[2] - 10)/5 <= 1)
    @addConstraint(m, -(D[1] - 20)/10 - (D[2] - 10)/5 <= 1)
    @test solve(m, prefer_cuts=cuts, add_box=cuts?1e2:false) == :Optimal
    @test isapprox(getValue(B[1]), 10.0, atol=TOL)
    @test isapprox(getValue(B[2]),  5.0, atol=TOL)

    # Add a deterministic constraint
    @addConstraint(m, B[1] <= 8)
    @test solve(m, prefer_cuts=cuts, add_box=cuts?1e2:false) == :Optimal
    @test isapprox(getValue(B[1]), 8.0, atol=TOL)
    @test isapprox(getValue(B[2]), 5.0, atol=TOL)

    # Add an uncertain constraint (and disambiguate objective)
    @addConstraint(m, B[1] + B[2] <= (D[1] + D[2])/2)
    @setObjective(m, Max, 3.1*S[2] + 3.0*S[1] - sum(B))
    @test solve(m, prefer_cuts=cuts, add_box=cuts?1e2:false) == :Optimal
    @test isapprox(getValue(B[1]), 5.0, atol=TOL)
    @test isapprox(getValue(B[2]), 5.0, atol=TOL)

    if solver in soc_solvers
        # Add a do-nothing ellipse constraint too
        @addConstraint(m, 1000 >= norm(D))
        @test solve(m, prefer_cuts=cuts, add_box=cuts?1e2:false) == :Optimal
        @test isapprox(getValue(B[1]), 5.0, atol=TOL)
        @test isapprox(getValue(B[2]), 5.0, atol=TOL)
    end
end  # "Resolving with..."


@testset "show_cuts" begin
    m = RobustModel()
    @defVar(m, 0 <= x <= 10)
    @defUnc(m, 0 <= u <= 10)
    @setObjective(m, Max, 10x)
    @addConstraint(m, u*x <= 7)
    @addConstraint(m, u <= 7)
    old_stdout = STDOUT
    rd, wr = redirect_stdout()
    solve(m, prefer_cuts=true, show_cuts=true)
    redirect_stdout(old_stdout)
    @test isapprox(getValue(x), 1.0, atol=TOL)
end  #show_cuts...""

end  # "GeneralOracle"
