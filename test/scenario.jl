using JuMPeR
using Base.Test

let
    rm = RobustModel()
    @defVar(rm, x >= 0)
    @defUnc(rm, 3 <= u <= 5)
    @defUnc(rm, 1 <= v <= 1)
    @setObjective(rm, Max, x)
    addConstraint(rm, v*x <= u)

    # Now we get tricky and add a scenario that is actually
    # outside the uncertainty set - but we don't check that
    addScenario(rm, [u => 1.5, v => 3.0])
    
    # This scenario shouldn't be added because it isn't
    # able to complete define a constraint
    addScenario(rm, [v => 9999.0])

    solveRobust(rm)
    @test_approx_eq getValue(x) 0.5
end

let
    # Test retrieval of active cuts
    rm = RobustModel()
    @defVar(rm, x >= 0)
    @defUnc(rm, 3 <= u <= 5)
    @defUnc(rm, 1 <= v <= 1)
    @setObjective(rm, Max, x)
    mycon = addConstraint(rm, v*x <= u)
    solveRobust(rm, prefer_cuts=true, active_cuts=true)
    ascen = getScenario(mycon)
    @test_approx_eq getUncValue(ascen, u) 3.0
    @test_approx_eq getUncValue(ascen, v) 1.0
end