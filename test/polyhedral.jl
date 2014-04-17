
function Test1(pref)
  println("  Test1")
  m = RobustModel(solver=solver)
  @defVar(m, x[1:2] >= 0)
  @defUnc(m, 0.3 <= u <= 0.5)
  setObjective(m, :Max, x[1] + x[2])
  addConstraint(m, u*x[1] + 1*x[2] <= 2.)
  addConstraint(m, 1*x[1] + 1*x[2] <= 6.)
  status = solveRobust(m, prefer_cuts=pref)
  @test_approx_eq getValue(x[1]) 4.0
  @test_approx_eq getValue(x[2]) 0.0
end

function Test2(pref)
  println("  Test2")
  m = RobustModel(solver=solver)
  @defVar(m, x[1:2] >= 0)
  @defUnc(m, 0.3 <= u1 <= 0.5)
  @defUnc(m, 0.0 <= u2 <= 2.0)
  setObjective(m, :Max, x[1] + x[2])
  addConstraint(m, u1*x[1] + 1*x[2] <= 2.)
  addConstraint(m, u2*x[1] + 1*x[2] <= 6.)
  status = solveRobust(m, prefer_cuts=pref)
  @test_approx_eq getValue(x[1]) (2.0+2.0/3.0)
  @test_approx_eq getValue(x[2]) (    2.0/3.0)
end

function Test2Flip(pref)
  println("  Test2Flip")
  m = RobustModel(solver=solver)
  @defVar(m, x[1:2] >= 0)
  @defUnc(m, 0.3 <= u1 <= 0.5)
  @defUnc(m, 0.0 <= u2 <= 2.0)
  setObjective(m, :Max, x[1] + x[2])
  addConstraint(m, -1*u1*x[1] + -1*x[2] >= -2.)
  addConstraint(m, -1*u2*x[1] + -1*x[2] >= -6.)
  status = solveRobust(m, prefer_cuts=pref)
  @test_approx_eq getValue(x[1]) (2.0+2.0/3.0)
  @test_approx_eq getValue(x[2]) (    2.0/3.0)
end

function Test2IP(pref)
  println("  Test2IP")
  m = RobustModel(solver=solver)
  @defVar(m, x[1:2] >= 0, Int)
  @defUnc(m, 0.3 <= u1 <= 0.5)
  @defUnc(m, 0.0 <= u2 <= 2.0)
  setObjective(m, :Max, 1.1*x[1] + x[2])
  addConstraint(m, u1*x[1] + 1*x[2] <= 2.)
  addConstraint(m, u2*x[1] + 1*x[2] <= 6.)
  status = solveRobust(m, prefer_cuts=true)
  @test_approx_eq getValue(x[1]) 3.0
  @test_approx_eq getValue(x[2]) 0.0
end

function Test3(pref)
  println("  Test3")
  m = RobustModel(solver=solver)
  @defVar(m, x[1:2] >= 0)
  @defUnc(m, 0.3 <= u1 <= 1.5)
  @defUnc(m, 0.5 <= u2 <= 1.5)
  setObjective(m, :Max, x[1] + x[2])
  addConstraint(m, u1*x[1] <= 3)
  addConstraint(m, u2*x[2] <= 1)
  addConstraint(m, (2.0*u1-2.0) + (4.0*u2-2.0) <= +1)
  addConstraint(m, (2.0*u1-2.0) + (4.0*u2-2.0) >= -1)
	status = solveRobust(m, prefer_cuts=pref)
  @test_approx_eq getValue(x[1]) (2.0)
  @test_approx_eq getValue(x[2]) (10.0/11.0)
end

function Test3IP(pref)
  println("  Test3IP")
  m = RobustModel(solver=solver)
  @defVar(m, x[1:2] >= 0, Int)
  @defUnc(m, 0.3 <= u1 <= 1.5)
  @defUnc(m, 0.5 <= u2 <= 1.5)
  setObjective(m, :Max, x[1] + x[2])
  addConstraint(m, u1*x[1] <= 3)
  addConstraint(m, u2*x[2] <= 1)
  addConstraint(m, (2.0*u1-2.0) + (4.0*u2-2.0) <= +1)
  addConstraint(m, (2.0*u1-2.0) + (4.0*u2-2.0) >= -1)
  status = solveRobust(m, prefer_cuts=pref)
  @test_approx_eq getValue(x[1]) 2.0
  @test_approx_eq getValue(x[2]) 0.0
end

function Test4(pref)
  println("  Test4")
  m = RobustModel(solver=solver)
  @defVar(m, x >= 0)
  @defUnc(m, 3.0 <= u <= 4.0);
  setObjective(m, :Max, 1.0x)
  addConstraint(m, 1.0*x <= 1.0*u)
  status = solveRobust(m, prefer_cuts=pref)
  @test_approx_eq getValue(x) 3.0
end

function Test5(pref)
  println("  Test5")
  m = RobustModel()
  @defUnc(m, u[1:5] >= 0)
  for ix = 1:5
      addConstraint(m, u[ix] <= float(ix))
  end
  addConstraint(m, sum(u[:]) <= 5.)
  @defVar(m, t)
  addConstraint(m, t >= sum(u[:]))
  @setObjective(m, Min, t)
  solveRobust(m, prefer_cuts=pref)
  @test_approx_eq getObjectiveValue(m) 5.0
end

function Test6(pref, variant)
  println("  Test6 variant $variant")
  rm = RobustModel()
  @defUnc(rm, u >=0)
  addConstraint(rm, u <=0)
  @defVar(rm, x >=0)
  @defVar(rm, shed >=0)
  @setObjective(rm, Min, x + shed)
  variant == 0 && addConstraint(rm, x - u + 3.46 - shed <= 0)
  variant == 1 && addConstraint(rm, x - u + 3.46 <= shed)
  variant == 2 && addConstraint(rm, x - u <= shed - 3.46)
  variant == 3 && addConstraint(rm, x <= u + shed - 3.46)
  variant == 4 && addConstraint(rm, 0 <= -x + u + shed - 3.46)
  variant == 5 && addConstraint(rm, 3.46 <= -x + u + shed)
  variant == 6 && addConstraint(rm, 3.46 - shed <= -x + u)
  variant == 7 && addConstraint(rm, 3.46 + x <= shed + u)
  solveRobust(rm, prefer_cuts=pref)
  @test_approx_eq getObjectiveValue(rm) 3.46
end

for pref in [true,false]
  println(" prefer_cuts:", pref)
  Test1(pref)
  Test2(pref)
  Test2Flip(pref)
  Test2IP(pref)
  Test3(pref)
  Test3IP(pref)
  Test4(pref)
  Test5(pref)
  for variant = 0:7
    Test6(pref,variant)
  end
end