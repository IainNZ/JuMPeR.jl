using JuMPeR
using FactCheck
import JuMP.REPLMode, JuMP.IJuliaMode

# Helper function to test IO methods work correctly
function io_test(mode, obj, exp_str; repl=:both)
    if mode == REPLMode
        repl != :show  && @fact sprint(print, obj) => exp_str
        repl != :print && @fact sprint(show,  obj) => exp_str
    else
        @fact sprint(writemime, "text/latex", obj) => "\$\$ $exp_str \$\$"
    end
end

lastc(rm)  = conToStr(JuMPeR.getRobust(rm).uncertainconstr[end])
lastuc(rm) = conToStr(JuMPeR.getRobust(rm).uncertaintyset[end])


facts("[macro] JuMPContainer{Uncertain}") do
    le, ge = JuMP.repl_leq, JuMP.repl_geq

    m = RobustModel()
    @defUnc(m,      bnd_free[2:5])
    @defUnc(m,      bnd_lowb[2:5] >= 2)
    @defUnc(m,      bnd_high[2:5] <= 5)
    @defUnc(m, 2 <= bnd_both[2:5] <= 5)
    @defUnc(m,      bnd_difflo[i=2:5] >= i)
    @defUnc(m,      bnd_diffup[i=2:5] <= i)
    @defUnc(m, i <= bnd_diffbo[i=2:5] <= 2i)
    @defUnc(m, i <= bnd_difflo_with_up[i=2:5] <= 5)
    @defUnc(m, 2 <= bnd_diffup_with_lo[i=2:5] <= i)

    io_test(REPLMode, bnd_free, "bnd_free[i] free for all i in {2,3,4,5}")
    io_test(REPLMode, bnd_lowb, "bnd_lowb[i] $ge 2 for all i in {2,3,4,5}")
    io_test(REPLMode, bnd_high, "bnd_high[i] $le 5 for all i in {2,3,4,5}")
    io_test(REPLMode, bnd_both, "2 $le bnd_both[i] $le 5 for all i in {2,3,4,5}")
    io_test(REPLMode, bnd_difflo, "bnd_difflo[i] $ge .. for all i in {2,3,4,5}")
    io_test(REPLMode, bnd_diffup, "bnd_diffup[i] $le .. for all i in {2,3,4,5}")
    io_test(REPLMode, bnd_diffbo, ".. $le bnd_diffbo[i] $le .. for all i in {2,3,4,5}")
    io_test(REPLMode, bnd_difflo_with_up, ".. $le bnd_difflo_with_up[i] $le 5 for all i in {2,3,4,5}")
    io_test(REPLMode, bnd_diffup_with_lo, "2 $le bnd_diffup_with_lo[i] $le .. for all i in {2,3,4,5}")
end


#=
@addConstraint(rm, u*x <= 5)
@test lastc(rm) == "u x <= 5"

@addConstraint(rm, sum{v[i]*y[i], i=1:5; i!=3} <= 9)
@test lastc(rm) == "v[1] y[1] + v[2] y[2] + v[4] y[4] + v[5] y[5] <= 9"

@addConstraint(rm, sum{v[i], i=1:5} == 1)
@test lastuc(rm) == "v[5] + v[4] + v[3] + v[2] + v[1] == 1"

@addConstraint(rm, u == sum{i*v[i], i=1:3})
@test lastuc(rm) == "-3 v[3] - 2 v[2] - v[1] + u == 0"

@addConstraint(rm, sum{i*(u+v[i])*(y[i]+x), i=1:2:5} <= 0)
@test lastc(rm) == "(u + v[1]) y[1] + (u + v[1]) x + (3 u + 3 v[3]) y[3] + (3 u + 3 v[3]) x + (5 u + 5 v[5]) y[5] + (5 u + 5 v[5]) x <= 0"

end
=#