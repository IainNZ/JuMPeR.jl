using JuMP, JuMPeR
using Base.Test

#----------------------------------------------------------------------

println("  Test 1")
m = RobustModel()
r = JuMPeR.getRobust(m)
@defUnc(m, 0 <= u[1:2,1:4] <= 1)
for i = 1:2
    for j = 2:4
        addConstraint(m, u[i,1] + u[i,j] >= 0)
    end
end
ret = JuMPeR.detect_components(r.numUncs, r.uncertaintyset)
@test ret[1] == [1,1,1,1,2,2,2,2]
@test ret[2] == [1,1,1,2,2,2]

@defVar(m, x[1:2] >= 0)
@setObjective(m, Min, sum(x))
setDefaultOracle!(m, JuMPeR.GeneralGraphOracle())
addConstraint(m, x[1] >= u[1,1])
addConstraint(m, x[2] >= u[1,4])
solveRobust(m, prefer_cuts=true)    # CUTS
@test_approx_eq getValue(x[1]) 1.0
@test_approx_eq getValue(x[2]) 1.0
solveRobust(m, prefer_cuts=false)   # REFORM
@test_approx_eq getValue(x[1]) 1.0
@test_approx_eq getValue(x[2]) 1.0

#----------------------------------------------------------------------

println("  Test 2")
m = RobustModel()
r = JuMPeR.getRobust(m)
@defUnc(m, 2 <= u[1:3] <= 2)
addConstraint(m, u[1] + u[2]        >= 2)
addConstraint(m,        u[2] + u[3] >= 2)
ret = JuMPeR.detect_components(r.numUncs, r.uncertaintyset)
@test ret[1] == [1,1,1]
@test ret[2] == [1,1]

@defVar(m, x[1:2] >= 0)
@setObjective(m, Max, sum(x))
setDefaultOracle!(m, JuMPeR.GeneralGraphOracle())
addConstraint(m, x[1] <= u[1])
addConstraint(m, x[2] <= u[3])
solveRobust(m, prefer_cuts=true)    # CUTS
@test_approx_eq getValue(x[1]) 2.0
@test_approx_eq getValue(x[2]) 2.0
solveRobust(m, prefer_cuts=false)   # REFORM
@test_approx_eq getValue(x[1]) 2.0
@test_approx_eq getValue(x[2]) 2.0

#----------------------------------------------------------------------

println("  Test 3")
m = RobustModel()
r = JuMPeR.getRobust(m)
@defUnc(m, u[1:4] >= 0)
addConstraint(m, u[1]        + u[3]        == 1)
addConstraint(m,        u[2]               == 1)
addConstraint(m,        u[2] + u[3]        == 1)
addConstraint(m,                      u[4] <= 5)
ret = JuMPeR.detect_components(r.numUncs, r.uncertaintyset)
@test ret[1] == [1,1,1,2]
@test ret[2] == [1,1,1,2]

@defVar(m, x[1:2] >= 0)
@setObjective(m, Min, sum(x))
setDefaultOracle!(m, JuMPeR.GeneralGraphOracle())
addConstraint(m, x[1] >= u[1])
addConstraint(m, x[2] >= u[4])
printRobust(m)
solveRobust(m, prefer_cuts=true)    # CUTS
@test_approx_eq getValue(x[1]) 1.0
@test_approx_eq getValue(x[2]) 5.0
solveRobust(m, prefer_cuts=false)   # REFORM
@test_approx_eq getValue(x[1]) 1.0
@test_approx_eq getValue(x[2]) 5.0