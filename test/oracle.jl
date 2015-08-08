#############################################################################
# JuMPeR
# Julia for Mathematical Programming - extension for Robust Optimization
# See http://github.com/IainNZ/JuMPeR.jl
#############################################################################
# test/oracle.jl
# Test oracle interface & tools
#############################################################################

using JuMP, JuMPeR
using FactCheck
using Compat

if !(:lp_solvers in names(Main))
    println("Loading solvers...")
    include(joinpath(Pkg.dir("JuMP"),"test","solvers.jl"))
end
lp_solvers = filter(s->(string(typeof(s))!="SCS.SCSSolver"), lp_solvers)

type IncompleteOracle <: AbstractOracle end

facts("[oracle] Check interface throws") do
    @fact_throws registerConstraint(IncompleteOracle(), RobustModel(), 1, nothing)
    @fact_throws setup(IncompleteOracle(), RobustModel(), nothing)
    @fact_throws generateReform(IncompleteOracle(), Model(), RobustModel(), Int[])
    @fact_throws generateCut(IncompleteOracle(), Model(), RobustModel(), Int[])
end

# build_cut_objective
# build_cut_objective_sparse
# build_certain_constraint
# is_constraint_violated
facts("[oracle] Test oracle utilities") do
    rm = RobustModel()
    @defVar(rm, x[1:4] >= 0)
    @defUnc(rm, u[1:5])
    @addConstraint(rm,  (3*u[1] + 2.0) * x[1] +
                        (  u[2] - 1.0) * x[2] +
                        (u[1] +  u[3]) * x[3] +
                        (u[3] +2*u[4]) * x[4] <= 5.0 + u[5])
    col_val = [2.0, 3.0, 4.0, 5.0]
    unc_con = JuMPeR.getRobust(rm).uncertainconstr[end]

    # -------------------
    sense, unc_coeffs, lhs_const = JuMPeR.build_cut_objective(rm, unc_con, col_val)
    @fact sense --> :Max
    @fact unc_coeffs[1] --> roughly(3.0*2.0+4.0, 1e-6)
    @fact unc_coeffs[2] --> roughly(3.0, 1e-6)
    @fact unc_coeffs[3] --> roughly(4.0+5.0, 1e-6)
    @fact unc_coeffs[4] --> roughly(2*5.0, 1e-6)
    @fact unc_coeffs[5] --> roughly(-1.0, 1e-6)
    @fact lhs_const     --> roughly(2.0*2.0-1.0*3.0, 1e-6)

    # -------------------
    sense, unc_coeffs, lhs_const = JuMPeR.build_cut_objective_sparse(unc_con, col_val)
    sort!(unc_coeffs)
    
    @fact sense --> :Max
    for i = 1:5
        @fact unc_coeffs[i][1] --> i
    end
    @fact unc_coeffs[1][2] --> roughly(3.0*2.0+4.0, 1e-6)
    @fact unc_coeffs[2][2] --> roughly(3.0, 1e-6)
    @fact unc_coeffs[3][2] --> roughly(4.0+5.0, 1e-6)
    @fact unc_coeffs[4][2] --> roughly(2*5.0, 1e-6)
    @fact unc_coeffs[5][2] --> roughly(-1.0, 1e-6)
    @fact lhs_const        --> roughly(2.0*2.0-1.0*3.0, 1e-6)

    # -------------------
    unc_val = [1.0, 2.0, 3.0, 4.0, 5.0]
    new_con = JuMPeR.build_certain_constraint(rm, unc_con, unc_val)
    if VERSION < v"0.4.0-"
        @fact conToStr(new_con) --> "5 x[1] + x[2] + 4 x[3] + 11 x[4] $(JuMP.repl_leq) 10"
    else
        @fact conToStr(new_con) --> "x[2] + 5 x[1] + 4 x[3] + 11 x[4] $(JuMP.repl_leq) 10"
    end

    # Bit of a hack to test build from JuMPDict
    inner_m = Model(solver=lp_solvers[1])
    @defVar(inner_m, i <= inner_u[i=1:5] <= i)
    @setObjective(inner_m, Max, sum(inner_u))
    solve(inner_m)
    new_con = JuMPeR.build_certain_constraint(rm, unc_con, getValue(inner_u))
    if VERSION < v"0.4.0-"
        @fact conToStr(new_con) --> "5 x[1] + x[2] + 4 x[3] + 11 x[4] $(JuMP.repl_leq) 10"
    else
        @fact conToStr(new_con) --> "x[2] + 5 x[1] + 4 x[3] + 11 x[4] $(JuMP.repl_leq) 10"
    end

    # -------------------
    lhs_val = dot([5,1,4,11],[2,3,4,5])
    @fact JuMPeR.check_cut_status(new_con, lhs_val, 1e-6) --> :Violate
    lhs_val = dot([5,1,4,11],[2,0,0,0])
    @fact JuMPeR.check_cut_status(new_con, lhs_val, 1e+6) --> :Active
    lhs_val = dot([5,1,4,11],[0,0,0,0])
    @fact JuMPeR.check_cut_status(new_con, lhs_val, 1e-6) --> :Slack
end