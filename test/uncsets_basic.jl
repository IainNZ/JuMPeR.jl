#-----------------------------------------------------------------------
# JuMPeR  --  JuMP Extension for Robust Optimization
# http://github.com/IainNZ/JuMPeR.jl
#-----------------------------------------------------------------------
# Copyright (c) 2016: Iain Dunning
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
#-----------------------------------------------------------------------
# test/uncsets_basic.jl
# Test BasicUncertaintySet for general polyhedral uncertainty sets.
#-----------------------------------------------------------------------

using JuMP, JuMPeR
using Test

TOL = 1e-4

if !(:lp_solvers in names(Main))
    printstyled("Loading solvers...\n", color = :magenta)
    include(joinpath(dirname(pathof(JuMP)),"..","test","solvers.jl"))
end
lp_solvers  = filter(s->(!occursin("SCSSolver", string(typeof(s)))), lp_solvers)
soc_solvers = filter(s->(!occursin("SCSSolver", string(typeof(s)))), soc_solvers)
solver_name(solver) = split(string(typeof(solver)),".")[2]

@testset "BasicUncertaintySet polyhedral" begin
printstyled("BasicUncertaintySet polyhedral...\n", color = :yellow)
printstyled("  LP tests...\n", color = :yellow)
@testset "LPs with $(solver_name(solver)), cuts=$cuts" for
                        solver in lp_solvers, cuts in [true,false]

    @testset "Test 1" begin
        m = RobustModel(solver=solver)
        @variable(m, x[1:2] >= 0)
        @uncertain(m, 0.3 <= u <= 0.5)
        @objective(m, Max, x[1] + x[2])
        @constraint(m, u*x[1] + 1*x[2] <= 2)
        @constraint(m, 1*x[1] + 1*x[2] <= 6)
        @test solve(m, prefer_cuts=cuts) == :Optimal
        @test isapprox(getvalue(x[1]), 4.0, atol=TOL)
        @test isapprox(getvalue(x[2]), 0.0, atol=TOL)
    end  # "Test 1"

    @testset "Test 2" begin
        m = RobustModel(solver=solver)
        @variable(m, 0 <= x[1:2] <= 9)
        @uncertain(m, 0.3 <= u1 <= 0.5)
        @uncertain(m, 0.0 <= u2 <= 2.0)
        @objective(m, Max, x[1] + x[2])
        @constraint(m, u1*x[1] + 1*x[2] <= 2)
        @constraint(m, u2*x[1] + 1*x[2] <= 6)
        @test solve(m, prefer_cuts=cuts) == :Optimal
        @test isapprox(getvalue(x[1]), 2.0+2.0/3.0, atol=TOL)
        @test isapprox(getvalue(x[2]),     2.0/3.0, atol=TOL)
    end  # "Test 2"

    @testset "Test 3" begin
        m = RobustModel(solver=solver)
        @variable(m, 0 <= x[1:2] <= 9)
        @uncertain(m, 0.3 <= u1 <= 1.5)
        @uncertain(m, 0.5 <= u2 <= 1.5)
        @objective(m, Max, x[1] + x[2])
        @constraint(m, u1*x[1] <= 3)
        @constraint(m, u2*x[2] <= 1)
        @constraint(m, (2.0*u1-2.0) + (4.0*u2-2.0) <= +1)
        @constraint(m, (2.0*u1-2.0) + (4.0*u2-2.0) >= -1)
        @test solve(m, prefer_cuts=cuts) == :Optimal
        @test isapprox(getvalue(x[1]), 2.0,       atol=TOL)
        @test isapprox(getvalue(x[2]), 10.0/11.0, atol=TOL)
    end  # "Test 3"

    @testset "Test 4" begin
        m = RobustModel(solver=solver)
        @variable(m, 0 <= x <= 9)
        @uncertain(m, u <= 4.0)
        @constraint(m, u >= 3.0)
        @objective(m, Max, 1.0x)
        @constraint(m, x <= u)
        @test solve(m, prefer_cuts=cuts) == :Optimal
        @test isapprox(getvalue(x), 3.0, atol=TOL)
    end  # "Test 4"

    @testset "Test 5" begin
        m = RobustModel(solver=solver)
        @uncertain(m, u[1:5] >= 0)
        for ix = 1:5
            @constraint(m, u[ix] <= float(ix))
        end
        @constraint(m, sum(u) <= 5)
        @variable(m, t >= 0)
        @constraint(m, t >= sum(u))
        @objective(m, Min, t)
        @test solve(m, prefer_cuts=cuts) == :Optimal
        @test isapprox(getobjectivevalue(m), 5.0, atol=TOL)
    end  # "Test 5"

    @testset "Test 6 variant $variant" for variant in 0:7
        rm = RobustModel(solver=solver)
        @uncertain(rm, u >= 0)
        @constraint(rm, u <= 0)
        @variable(rm, x >= 0)
        @variable(rm, shed >= 0)
        @objective(rm, Min, x + shed)
        variant == 0 && @constraint(rm, x - u + 3.46 - shed <= 0)
        variant == 1 && @constraint(rm, x - u + 3.46 <= shed)
        variant == 2 && @constraint(rm, x - u <= shed - 3.46)
        variant == 3 && @constraint(rm, x <= u + shed - 3.46)
        variant == 4 && @constraint(rm, 0 <= -x + u + shed - 3.46)
        variant == 5 && @constraint(rm, 3.46 <= -x + u + shed)
        variant == 6 && @constraint(rm, 3.46 - shed <= -x + u)
        variant == 7 && @constraint(rm, 3.46 + x <= shed + u)
        @test solve(rm, prefer_cuts=cuts) == :Optimal
        @test isapprox(getobjectivevalue(rm), 3.46, atol=TOL)
    end  # "Test 6"

    @testset "Test 7" begin
        m = RobustModel(solver=solver)
        @variable(m, 0 <= x <= 9)
        @uncertain(m, 0.5 <= u <= 0.5)
        @objective(m, Max, x)
        @constraint(m, u*x + u <= 2)
        @test solve(m, prefer_cuts=cuts) == :Optimal
        @test isapprox(getvalue(x), 3.0, atol=TOL)
    end  # "Test 7"

    @testset "Test 8 (unbounded LP)" begin
        m = RobustModel(solver=solver)
        @variable(m, x >= 0)
        @uncertain(m, 8 <= u <= 9)
        @objective(m, Max, x)
        @constraint(m, u*x >= 25)
        @test solve(m, prefer_cuts=cuts, suppress_warnings=true) == :Unbounded
    end  # "Test 8"

    if !occursin("IpoptSolver", "$(typeof(solver))")  # reports UserLimit
    @testset "Test 9 (infeasible LP)" begin
        m = RobustModel(solver=solver)
        @variable(m, x)
        @uncertain(m, 8 <= u <= 9)
        @objective(m, Max, x)
        @constraint(m, x == 3)
        @constraint(m, u*x <= 25)
        @test solve(m, prefer_cuts=cuts, suppress_warnings=true) == :Infeasible
    end; end  # "Test 9"

    @testset "Test 10 (empty unc. set)" begin
        m = RobustModel(solver=solver)
        @variable(m, x >= 0)
        @uncertain(m, 8 <= u <= 7)
        @objective(m, Min, x)
        @constraint(m, u*x >= 25)
        if cuts
            # Cutting plane problem is infeasible, so an error will be
            # thrown that we catch here.
            @test_throws Exception solve(m, prefer_cuts=cuts)
        else
            # Reformulation doesn't care, although it is a bit weird
            @test solve(m, prefer_cuts=cuts) == :Optimal
        end
    end

    if !occursin("IpoptSolver", "$(typeof(solver))")  # reports UserLimit
    @testset "Test 11 (unbounded unc. set)" begin
        m = RobustModel(solver=solver)
        @variable(m, x >= 1)
        @uncertain(m, u <= 7)
        @objective(m, Min, x)
        @constraint(m, u*x >= 25)
        if cuts
            # Cutting plane problem is unbounded, so an error will be
            # thrown that we catch here
            @test_throws Exception solve(m, prefer_cuts=cuts)
        else
            @test solve(m, prefer_cuts=false, suppress_warnings=true) == :Infeasible
        end
    end; end

end  # "LPs with ..."

printstyled("  MILP tests...\n", color = :yellow)
@testset "MILPs with $(solver_name(solver)), cuts=$cuts" for
                        solver in lazy_solvers, cuts in [true,false]

    @testset "Test 1" begin
        m = RobustModel(solver=solver)
        @variable(m, 0 <= x[1:2] <= 9, Int)
        @uncertain(m, 0.3 <= u1 <= 0.5)
        @uncertain(m, 0.0 <= u2 <= 2.0)
        @objective(m, Max, 1.1*x[1] + x[2])
        @constraint(m, u1*x[1] + 1*x[2] <= 2)
        @constraint(m, u2*x[1] + 1*x[2] <= 6)
        @test solve(m, prefer_cuts=cuts) == :Optimal
        @test isapprox(getvalue(x[1]), 3.0, atol=TOL)
        # @test isapprox(getvalue(x[2]), 0.0, atol=TOL) # flakey test due to bug in GLPK
    end

    @testset "Test 2" begin
        m = RobustModel(solver=solver)
        @variable(m, 0 <= x[1:2] <= 9, Int)
        @uncertain(m, 0.3 <= u1 <= 1.5)
        @uncertain(m, 0.5 <= u2 <= 1.5)
        @objective(m, Max, x[1] + x[2])
        @constraint(m, u1*x[1] <= 3)
        @constraint(m, u2*x[2] <= 1)
        @constraint(m, (2*u1-2) + (4*u2-2) <= +1)
        @constraint(m, (2*u1-2) + (4*u2-2) >= -1)
        @test solve(m, prefer_cuts=cuts) == :Optimal
        @test isapprox(getvalue(x[1]), 2.0, atol=TOL)
        @test isapprox(getvalue(x[2]), 0.0, atol=TOL)
    end

    if cuts
    @testset "Test 3 (integer unc. set)" begin
        rm = RobustModel(solver=solver)
        @variable(rm, 0 <= x <= 2)
        @uncertain(rm, 0.5 <= u <= 1.5, Int)
        @objective(rm, Max, x)
        @constraint(rm, x <= u)
        @test solve(rm, prefer_cuts=true) == :Optimal
        @test isapprox(getvalue(x), 1.0, atol=TOL)
    end; end

    @testset "Test 4 (infeasible MILP)" begin
        m = RobustModel(solver=solver)
        @variable(m, x, Int)
        @uncertain(m, 8 <= u <= 9)
        @objective(m, Max, x)
        @constraint(m, x == 3)
        @constraint(m, u*x <= 25)  # u = 9, 9*3 = 27 > 25
        @test solve(m, prefer_cuts=cuts, suppress_warnings=true) == :Infeasible
    end

    if !occursin("GLPK", "$(typeof(solver))")  # reports Error
    @testset "Test 5 (unbounded MILP)" begin
        m = RobustModel(solver=solver)
        @variable(m, x >= 0, Int)
        @uncertain(m, 8 <= u <= 9)
        @objective(m, Max, x)
        @constraint(m, u*x >= 25)  # u = 8, x >= 0, > 25
        @test solve(m, prefer_cuts=cuts, suppress_warnings=true) == :Unbounded
    end; end

end  # "MILPs with..."


@testset "Scenarios with $(solver_name(solver)), cuts=$cuts" for
                        solver in lp_solvers, cuts in [true,false]
    @testset "Test 1" begin
        m = RobustModel(solver=solver)
        @variable(m, 0 <= x[1:2] <= 9)
        @uncertain(m, 0.3 <= u1 <= 0.5)
        @uncertain(m, 0.0 <= u2 <= 2.0)
        @objective(m, Max, x[1] + x[2])
        con1 = @constraint(m, u1*x[1] + 1*x[2] <= 2)
        con2 = @constraint(m, u2*x[1] + 1*x[2] <= 6)
        @test solve(m, prefer_cuts=cuts, active_scenarios=true) == :Optimal
        @test isapprox(getvalue(x[1]), 2.0+2.0/3.0, atol=TOL)
        @test isapprox(getvalue(x[2]),     2.0/3.0, atol=TOL)
        scen1 = getscenario(con1)
        @test isapprox(uncvalue(scen1, u1), 0.5, atol=TOL)
        scen2 = getscenario(con2)
        @test isapprox(uncvalue(scen2, u2), 2.0, atol=TOL)
    end
end


end  # "BasicUncertaintySet"
