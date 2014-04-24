using JuMPeR
using Base.Test

# build_certain_constraint
let
    rm = RobustModel()
    @defVar(rm, x[1:5] >= 0)
    @defUnc(rm, u[1:5])
    unc_con = (3*u[1] + 2.0) * x[1] +
              (  u[2] - 1.0) * x[2] +
                               x[3] +
              (u[3] + 2*u[4])* x[4] <=
              5.0 + u[5]
    unc_val = [1.0, 2.0, 3.0, 4.0, 5.0]
    new_con = JuMPeR.build_certain_constraint(unc_con, unc_val)
    @test conToStr(new_con) == "5 x[1] + x[2] + x[3] + 11 x[4] <= 10"
end
