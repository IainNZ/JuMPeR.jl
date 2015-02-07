using JuMP, JuMPeR
using Base.Test

type IncompleteOracle <: AbstractOracle end

# build_cut_objective
# build_cut_objective_sparse
# build_certain_constraint
# is_constraint_violated
let

    @test_throws ErrorException registerConstraint(IncompleteOracle(), RobustModel(), 1, nothing)
    @test_throws ErrorException setup(IncompleteOracle(), RobustModel(), nothing)
    @test_throws ErrorException generateReform(IncompleteOracle(), Model(), RobustModel(), Int[])
    @test_throws ErrorException generateCut(IncompleteOracle(), Model(), RobustModel(), Int[])


    # -------------------
    rm = RobustModel()
    @defVar(rm, x[1:4] >= 0)
    @defUnc(rm, u[1:5])
    unc_con = (3*u[1] + 2.0) * x[1] +
              (  u[2] - 1.0) * x[2] +
              (u[1] +  u[3]) * x[3] +
              (u[3] +2*u[4]) * x[4] <=
              5.0 + u[5]
    col_val = [2.0, 3.0, 4.0, 5.0]

    # -------------------
    sense, unc_coeffs, lhs_const = JuMPeR.build_cut_objective(rm, unc_con, col_val)
    @test sense == :Max
    @test_approx_eq unc_coeffs[1] 3.0*2.0+4.0
    @test_approx_eq unc_coeffs[2] 3.0
    @test_approx_eq unc_coeffs[3] 4.0+5.0
    @test_approx_eq unc_coeffs[4] 2*5.0
    @test_approx_eq unc_coeffs[5] -1.0
    @test_approx_eq lhs_const     2.0*2.0-1.0*3.0

    # -------------------
    sense, unc_coeffs, lhs_const = JuMPeR.build_cut_objective_sparse(unc_con, col_val)
    sort!(unc_coeffs)
    
    @test sense == :Max
    for i = 1:5
        @test unc_coeffs[i][1] == i
    end
    @test_approx_eq  unc_coeffs[1][2]  3.0*2.0+4.0
    @test_approx_eq  unc_coeffs[2][2]  3.0
    @test_approx_eq  unc_coeffs[3][2]  4.0+5.0
    @test_approx_eq  unc_coeffs[4][2]  2*5.0
    @test_approx_eq  unc_coeffs[5][2]  -1.0
    @test_approx_eq  lhs_const         2.0*2.0-1.0*3.0

    # -------------------
    unc_val = [1.0, 2.0, 3.0, 4.0, 5.0]
    new_con = JuMPeR.build_certain_constraint(rm, unc_con, unc_val)
    @test conToStr(new_con) == "5 x[1] + x[2] + 4 x[3] + 11 x[4] $(JuMP.repl_leq) 10"

    # Bit of a hack to test build from JuMPDict
    inner_m = Model()
    @defVar(inner_m, i <= inner_u[i=1:5] <= i)
    @setObjective(inner_m, Max, sum(inner_u))
    solve(inner_m)
    new_con = JuMPeR.build_certain_constraint(rm, unc_con, getValue(inner_u))
    @test conToStr(new_con) == "5 x[1] + x[2] + 4 x[3] + 11 x[4] $(JuMP.repl_leq) 10"

    # -------------------
    lhs_val = dot([5,1,4,11],[2,3,4,5])
    @test JuMPeR.check_cut_status(new_con, lhs_val, 1e-6) == :Violate
    lhs_val = dot([5,1,4,11],[2,0,0,0])
    @test JuMPeR.check_cut_status(new_con, lhs_val, 1e+6) == :Active
    lhs_val = dot([5,1,4,11],[0,0,0,0])
    @test JuMPeR.check_cut_status(new_con, lhs_val, 1e-6) == :Slack
end