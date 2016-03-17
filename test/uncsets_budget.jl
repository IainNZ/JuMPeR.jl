#-----------------------------------------------------------------------
# JuMPeR  --  JuMP Extension for Robust Optimization
# http://github.com/IainNZ/JuMPeR.jl
#-----------------------------------------------------------------------
# Copyright (c) 2016: Iain Dunning
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
#-----------------------------------------------------------------------
# test/uncsets_budget.jl
# Test the BudgetUncertaintySet
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

@testset "BudgetUncertaintySet" begin
print_with_color(:yellow, "BudgetUncertaintySet...\n")
@testset "LPs with $(solver_name(solver))" for solver in lp_solvers

    @testset "Test 1 (+x, +coeff)" begin
        μ, σ = [2.0, 6.0], [1.0, 1.0]
        m = RobustModel(solver=solver)
        @defVar(m, 0 <= x[1:2] <= 10)
        @defUnc(m, u[i=1:2])
        uncertainty_set!(m, JuMPeR.BudgetUncertaintySet(1, μ, σ))
        @setObjective(m, Max, 0.1*x[1] + 9.9*x[2])
        @addConstraint(m, u[1]*x[1] + u[2]*x[2] <= 21)
        @test solve(m) == :Optimal
        @test isapprox(getValue(x[1]), 0.0, atol=TOL)
        @test isapprox(getValue(x[2]), 3.0, atol=TOL)
    end

    @testset "Test 2 (+x, -coeff)" begin
        μ, σ = -[2.0, 6.0], [1.0, 1.0]
        m = RobustModel(solver=solver)
        @defVar(m, 0 <= x[1:2] <= 10)
        @defUnc(m, u[i=1:2])
        uncertainty_set!(m, JuMPeR.BudgetUncertaintySet(1, μ, σ))
        @setObjective(m, Max, 0.1*x[1] + 9.9*x[2])
        @addConstraint(m, u[1]*x[1] + u[2]*x[2] >= -21)
        @test solve(m) == :Optimal
        @test isapprox(getValue(x[1]), 0.0, atol=TOL)
        @test isapprox(getValue(x[2]), 3.0, atol=TOL)
    end

    @testset "Test 3 (-x, +coeff)" begin
        μ, σ = [2.0, 6.0], [1.0, 1.0]
        m = RobustModel(solver=solver)
        @defVar(m, -10 <= x[1:2] <= 0)
        @defUnc(m, u[i=1:2])
        uncertainty_set!(m, JuMPeR.BudgetUncertaintySet(1, μ, σ))
        @setObjective(m, Max, -0.1*x[1] - 9.9*x[2])
        @addConstraint(m, u[1]*x[1] + u[2]*x[2] >= -21)
        @test solve(m) == :Optimal
        @test isapprox(getValue(x[1]), 0.0, atol=TOL)
        @test isapprox(getValue(x[2]),-3.0, atol=TOL)
    end

    @testset "Test 4 (-x, -coeff)" begin
        μ, σ = -[2.0, 6.0], [1.0, 1.0]
        m = RobustModel(solver=solver)
        @defVar(m, -10 <= x[1:2] <= 0)
        @defUnc(m, u[i=1:2])
        uncertainty_set!(m, JuMPeR.BudgetUncertaintySet(1, μ, σ))
        @setObjective(m, Max, -0.1*x[1] - 9.9*x[2])
        @addConstraint(m, u[1]*x[1] + u[2]*x[2] <= 21)
        @test solve(m) == :Optimal
        @test isapprox(getValue(x[1]), 0.0, atol=TOL)
        @test isapprox(getValue(x[2]),-3.0, atol=TOL)
    end  # "Test 4"

    @testset "Test 5 (multi-variable u, positive)" begin
        m = RobustModel(solver=solver)
        @defVar(m, x ==  1)
        @defVar(m, y ==  2)
        @defVar(m, z == -3)
        @defVar(m, slack <= 100)
        @defUnc(m, u)
        uncertainty_set!(m, JuMPeR.BudgetUncertaintySet(1, [0.0], [3.0]))
        @setObjective(m, Max, slack)
        # (2y - 3z)u = (2*2 - 3*-3)u = (4 + 9)u = 13u
        # slack = 50 - x - 13u = 49-13u = 49-39 = 10
        @addConstraint(m, x + 2*u*y - 3*u*z + slack <= 50)
        @test solve(m) == :Optimal
        @test isapprox(getValue(slack), 10.0, atol=TOL)
    end

    @testset "Test 6 (multi-variable u, negative)" begin
        m = RobustModel(solver=solver)
        @defVar(m, x == -1)
        @defVar(m, y == -2)
        @defVar(m, z == +3)
        @defVar(m, slack <= 100)
        @defUnc(m, u)
        uncertainty_set!(m, JuMPeR.BudgetUncertaintySet(1, [0.0], [3.0]))
        @setObjective(m, Max, slack)
        # (2y - 3z)u = (2*-2 - 3*+3)u = (-4 - 9)u = -13u
        # slack = 50 - x + 13u = 51+13u = 51-39u = 12
        @addConstraint(m, x + 2*u*y - 3*u*z + slack <= 50)
        @test solve(m) == :Optimal
        @test isapprox(getValue(slack), 12.0, atol=TOL)
    end

    @testset "Test 7 (set mixing default), cuts=$cuts" for cuts in [true,false]
        m = RobustModel(solver=solver)
        @defVar(m, x <= 10)
        @defVar(m, y <= 10)
        @defUnc(m, u)
        @setObjective(m, Max, x + y)
        # Uses default
        @addConstraint(m, u <= 10)
        @addConstraint(m, u*x <= 50)
        # Use Budget
        us = JuMPeR.BudgetUncertaintySet(1,[5.],[5.])
        @addConstraint(m, u*y <= 50, uncset=us)
        @test solve(m, prefer_cuts=cuts) == :Optimal
        @test isapprox(getValue(x), 5.0, atol=TOL)
        @test isapprox(getValue(y), 5.0, atol=TOL)
    end

    @testset "Test 8 (set mixing explicit), cuts=$cuts" for cuts in [true,false]
        m = RobustModel(solver=solver)
        @defVar(m, x <= 10)
        @defVar(m, y <= 10)
        @defUnc(m, u)
        @setObjective(m, Max, x + y)
        # Use Basic
        us = JuMPeR.BasicUncertaintySet()
        @addConstraint(us, u <= 10)
        @addConstraint(m, u*x <= 50, uncset=us)
        # Use Budget
        us = JuMPeR.BudgetUncertaintySet(1,[5.],[5.])
        @addConstraint(m, u*y <= 50, uncset=us)
        @test solve(m, prefer_cuts=cuts) == :Optimal
        @test isapprox(getValue(x), 5.0, atol=TOL)
        @test isapprox(getValue(y), 5.0, atol=TOL)
    end

    @testset "Test 8 (set mixing, scens), cuts=$cuts" for cuts in [true,false]
        m = RobustModel(solver=solver)
        @defVar(m, x <= 10)
        @defVar(m, y <= 10)
        @defVar(m, z <= 10)
        @defUnc(m, u)
        @setObjective(m, Max, x + y + z)
        us1 = JuMPeR.BasicUncertaintySet();  @addConstraint(us1, u <= 10)
        cn1 = @addConstraint(m, u*x <= 50, uncset=us1)
        us2 = JuMPeR.BudgetUncertaintySet(1,[50.],[50.])
        cn2 = @addConstraint(m, u*y <= 500, uncset=us2)
        us3 = JuMPeR.BasicUncertaintySet();  @addConstraint(us3, u <= 1000)
        cn3 = @addConstraint(m, u*z <= 5000, uncset=us3)
        @test solve(m, prefer_cuts=cuts, active_cuts=true) == :Optimal
        @test isapprox(getValue(x), 5.0, atol=TOL)
        @test isapprox(getValue(y), 5.0, atol=TOL)
        @test isapprox(getValue(z), 5.0, atol=TOL)
        @test isapprox(unc_value(get(getScenario(cn1)), u), 10.0, atol=TOL)
        @test isapprox(unc_value(get(getScenario(cn2)), u), 100.0, atol=TOL)
        @test isapprox(unc_value(get(getScenario(cn3)), u), 1000.0, atol=TOL)
    end
end  # "LPs with"
end  # "BudgetUncertaintySet"
