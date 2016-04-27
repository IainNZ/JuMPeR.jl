#-----------------------------------------------------------------------
# JuMPeR  --  JuMP Extension for Robust Optimization
# http://github.com/IainNZ/JuMPeR.jl
#-----------------------------------------------------------------------
# Copyright (c) 2016: Iain Dunning
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
#-----------------------------------------------------------------------
# test/uncsets.jl
# Test UncertainySet interface & tools.
#-----------------------------------------------------------------------

using JuMP, JuMPeR
using BaseTestNext

if !(:lp_solvers in names(Main))
    print_with_color(:magenta, "Loading solvers...\n")
    include(joinpath(Pkg.dir("JuMP"),"test","solvers.jl"))
end
lp_solvers = filter(s->(!contains(string(typeof(s)),"SCSSolver")), lp_solvers)

@testset "UncertainySet" begin
print_with_color(:yellow, "UncertainySet...\n")

@testset "Check interface throws" begin
    eval(:(type IncompleteSet <: JuMPeR.AbstractUncertaintySet end))  # In global scope
    @test_throws ErrorException JuMPeR.setup_set(IncompleteSet(), RobustModel(), Int[], false, nothing)
    @test_throws ErrorException JuMPeR.generate_reform(IncompleteSet(), RobustModel(), Int[])
    @test_throws ErrorException JuMPeR.generate_cut(IncompleteSet(), RobustModel(), Int[])
    @test_throws ErrorException JuMPeR.generate_scenario(IncompleteSet(), RobustModel(), Int[])
end

# build_cut_objective
# build_cut_objective_sparse
# build_certain_constraint
# is_constraint_violated
@testset "Utilities" begin
    rm = RobustModel()
    @variable(rm, x[1:4] >= 0)
    @uncertain(rm, u[1:5])
    @constraint(rm,  (3*u[1] + 2.0) * x[1] +   # u1: 3*x1 = 3*2 = 6, c: 2*x1 = 2*2 =4
                        (  u[2] - 1.0) * x[2] +   # u2: 1*x2 = 1*3 = 3, c: -1*x2 = -1*3 = -3
                        (u[1] +  u[3]) * x[3] +   # u1: 3*x3 = 3*4 = 12, u3: 1*x3 = 4
                        (u[3] +2*u[4]) * x[4] <=  # u3: 1*x4 = 1*5 = 5, u4: 2*x4 = 10
                        5.0 + u[5])  # u5
    col_val = [2.0, 3.0, 4.0, 5.0]
    unc_con = JuMPeR.get_robust(rm).unc_constraints[end]

    # Accumulate the coefficients for each uncertain parameter using
    # the values of the decision variables and the coefficients, to build
    # the objective function for a cut generating problem
    sense, unc_coeffs, lhs_const = JuMPeR.build_cut_objective(rm, unc_con, col_val)
    @test sense == :Max
    @test unc_coeffs[1] == 6 + 4
    @test unc_coeffs[2] == 3
    @test unc_coeffs[3] == 4 + 5
    @test unc_coeffs[4] == 10
    @test unc_coeffs[5] == -1
    @test lhs_const     == 4 - 3

    # Same as above, but "sparse"
    sense, unc_coeffs, lhs_const = JuMPeR.build_cut_objective_sparse(unc_con, col_val)
    sort!(unc_coeffs)
    @test sense == :Max
    for i in 1:5
        @test unc_coeffs[i][1] == i
    end
    @test unc_coeffs[1][2] == 6 + 4
    @test unc_coeffs[2][2] == 3
    @test unc_coeffs[3][2] == 4 + 5
    @test unc_coeffs[4][2] == 10
    @test unc_coeffs[5][2] == -1
    @test lhs_const        == 4 - 3

    # -------------------
    unc_val = [1.0, 2.0, 3.0, 4.0, 5.0]
    new_con = JuMPeR.build_certain_constraint(unc_con, unc_val)
    @test string(new_con) == "5 x[1] + x[2] + 4 x[3] + 11 x[4] $(JuMP.repl[:leq]) 10"

    # Bit of a hack to test build from JuMP.JuMPDict
    inner_m = Model(solver=lp_solvers[1])
    @variable(inner_m, i <= inner_u[i=1:5] <= i)
    @objective(inner_m, Max, sum(inner_u))
    solve(inner_m)
    new_con = JuMPeR.build_certain_constraint(unc_con, getvalue(inner_u))
    @test string(new_con) == "5 x[1] + x[2] + 4 x[3] + 11 x[4] $(JuMP.repl[:leq]) 10"

    # -------------------
    lhs_val = 1.0*dot([5,1,4,11],[2,3,4,5])
    @test JuMPeR.check_cut_status(new_con, lhs_val, 1e-6) == :Violate
    lhs_val = 1.0*dot([5,1,4,11],[2,0,0,0])
    @test JuMPeR.check_cut_status(new_con, lhs_val, 1e+6) == :Active
    lhs_val = 1.0*dot([5,1,4,11],[0,0,0,0])
    @test JuMPeR.check_cut_status(new_con, lhs_val, 1e-6) == :Slack
end  # "Utilities"

end  # "UncertainySet"
