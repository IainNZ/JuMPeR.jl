#############################################################################
# JuMPeR
# Julia for Mathematical Programming - extension for Robust Optimization
# See http://github.com/IainNZ/JuMPeR.jl
#############################################################################
# test/print.jl
# Testing for all pretty-printing-related functionality
#############################################################################

using JuMP, JuMPeR
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

facts("[print] JuMPContainer{Uncertain}") do
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
