BertSimOracle()

function BSTest1(cuts)
    println(" BSTest1 $cuts - postive x, positive coeff")
    # Data
    n           = 2
    weight_low  = [1.0, 5.0]
    weight_high = [3.0, 7.0]
    values      = [0.1, 9.9]

    m = RobustModel(solver=solver)
    setDefaultOracle!(m, BertSimOracle(1))
    
    @defVar(m, x[1:n] >= 0)
    @defUnc(m, weight_low[i] <= u[i=1:n] <= weight_high[i])

    @setObjective(m, Max, sum{ values[i] * x[i], i=1:n})

    addConstraint(m, sum([u[i]*x[i] for i=1:n]) <= 21)
    
    status = solveRobust(m, prefer_cuts=cuts)
    @test_approx_eq getValue(x[1]) 0.0
    @test_approx_eq getValue(x[2]) 3.0
end

function BSTest2(cuts)
    println(" BSTest2 $cuts - postive x, negative coeff")
    # Data
    n           = 2
    weight_high = [-1.0, -5.0]
    weight_low  = [-3.0, -7.0]
    values      = [0.1, 9.9]

    m = RobustModel(solver=solver)
    
    @defVar(m, x[1:n] >= 0)
    @defUnc(m, weight_low[i] <= u[i=1:n] <= weight_high[i])

    @setObjective(m, Max, sum{ values[i] * x[i], i=1:n})

    addConstraint(m, sum([u[i]*x[i] for i=1:n]) >= -21, BertSimOracle(1))
    
    status = solveRobust(m, prefer_cuts=cuts)
    @test_approx_eq getValue(x[1]) 0.0
    @test_approx_eq getValue(x[2]) 3.0
end

function BSTest3(cuts)
    println(" BSTest3 $cuts - negative x, positive coeff")
    # Data
    n           = 2
    weight_low  = [1.0, 5.0]
    weight_high = [3.0, 7.0]
    values      = [-0.1, -9.9]

    m = RobustModel(solver=solver)
    
    @defVar(m, x[1:n] <= 0)
    @defUnc(m, weight_low[i] <= u[i=1:n] <= weight_high[i])

    @setObjective(m, Max, sum{ values[i] * x[i], i=1:n})

    addConstraint(m, sum([u[i]*x[i] for i=1:n]) >= -21, BertSimOracle(1))
    
    status = solveRobust(m, prefer_cuts=cuts)
    @test_approx_eq getValue(x[1]) 0.0
    @test_approx_eq getValue(x[2]) -3.0
end


function BSTest4(cuts)
    println(" BSTest4 $cuts - negative x, negative coeff")
    # Data
    n           = 2
    weight_high = [-1.0, -5.0]
    weight_low  = [-3.0, -7.0]
    values      = [-0.1, -9.9]

    m = RobustModel(solver=solver)
    
    @defVar(m, x[1:n] <= 0)
    @defUnc(m, weight_low[i] <= u[i=1:n] <= weight_high[i])

    @setObjective(m, Max, sum{ values[i] * x[i], i=1:n})

    addConstraint(m, sum([u[i]*x[i] for i=1:n]) <= 21, BertSimOracle(1))
    
    status = solveRobust(m, prefer_cuts=cuts)
    @test_approx_eq getValue(x[1]) 0.0
    @test_approx_eq getValue(x[2]) -3.0
end


BSTest1(true)
BSTest2(true)
BSTest3(true)
BSTest4(true)