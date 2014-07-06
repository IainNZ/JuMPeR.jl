function resolve_1(pref)
    m = RobustModel(solver=solver)
    @defVar(m, B[1:2] >= 0)
    @defVar(m, S[1:2] >= 0)
    
    # Uncertainty set
    @defUnc(m, D[1:2])
    addConstraint(m,  (D[1] - 20)/10 + (D[2] - 10)/5 <= 1.5)
    addConstraint(m,  (D[1] - 20)/10 - (D[2] - 10)/5 <= 1.5)
    addConstraint(m, -(D[1] - 20)/10 + (D[2] - 10)/5 <= 1.5)
    addConstraint(m, -(D[1] - 20)/10 - (D[2] - 10)/5 <= 1.5)

    # Constraints
    for i in 1:2
        addConstraint(m, S[i] <= D[i])
        addConstraint(m, S[i] <= B[i])
    end

    # Objective
    setObjective(m, :Max, 3*sum(S) - sum(B))

    # First solve
    solveRobust(m, prefer_cuts=pref)
    @test_approx_eq getValue(B[1]) 5.0
    @test_approx_eq getValue(B[2]) 2.5

    # Solve again with no changes
    solveRobust(m, prefer_cuts=pref)
    @test_approx_eq getValue(B[1]) 5.0
    @test_approx_eq getValue(B[2]) 2.5

    # Tighten uncertainty set 
    addConstraint(m,  (D[1] - 20)/10 + (D[2] - 10)/5 <= 1.0)
    addConstraint(m,  (D[1] - 20)/10 - (D[2] - 10)/5 <= 1.0)
    addConstraint(m, -(D[1] - 20)/10 + (D[2] - 10)/5 <= 1.0)
    addConstraint(m, -(D[1] - 20)/10 - (D[2] - 10)/5 <= 1.0)
    solveRobust(m, prefer_cuts=pref)
    @test_approx_eq getValue(B[1]) 10.0
    @test_approx_eq getValue(B[2])  5.0

    # Add a certain constraint
    addConstraint(m, B[1] <= 8)
    solveRobust(m, prefer_cuts=pref)
    @test_approx_eq getValue(B[1]) 8.0
    @test_approx_eq getValue(B[2]) 5.0

    # Add an uncertain constraint (and disambiguate objective)
    addConstraint(m, B[1] + B[2] <= (D[1] + D[2])/2)
    setObjective(m, :Max, 3.1*S[2] + 3.0*S[1] - sum(B))
    solveRobust(m, prefer_cuts=pref)
    @test_approx_eq getValue(B[1]) 5.0
    @test_approx_eq getValue(B[2]) 5.0
end

for cuts in [true, false]
    println(" prefer_cuts: $cuts")
    resolve_1(cuts)
end