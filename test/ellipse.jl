using JuMPeR
using Base.Test

# Check the construction of the constraint data
let
    rm = RobustModel()
    @defUnc(rm, u[1:5])

    vec = [1.0*u[1]+2.0*u[2]+1.0, 5.0*u[5]+1.0*u[1]+5.0, 4.0*u[4]+4.0]
    con = build_ellipse_constraint(vec, 2.0)
    @test con.F == [1.0 2.0 0.0 0.0;
                    1.0 0.0 5.0 0.0;
                    0.0 0.0 0.0 4.0]
    @test con.u == [  1,  2,  5,  4]
    @test con.g == [1.0, 5.0, 4.0]
    @test con.Gamma == 2.0

    vec = [u[5], u[3], u[1] + 9.0]
    con = build_ellipse_constraint(vec, 1.0)
    @test con.F == [1.0 0.0 0.0;
                    0.0 1.0 0.0;
                    0.0 0.0 1.0]
    @test con.u == [  5,  3,  1]
    @test con.g == [ 0., 0., 9.]
    @test con.Gamma == 1.0

    vec = [1.0*u[1]-5]
    con = build_ellipse_constraint(vec, 3.0)
    @test con.F == ones(1,1)
    @test con.u == [  1]
    @test con.g == [-5.]
    @test con.Gamma == 3.0

    vec = [2.0, 0.0]
    @test_throws ErrorException build_ellipse_constraint(vec, 0.0)

    addEllipseConstraint(rm, [u[5], u[3], u[1] + 9.0], 1.0)
end

# Use the ellipse constraint
if Pkg.installed("Gurobi") == nothing
    println("Cannot run ellipse uncertainty set solution tests - no Gurobi!")    
else
    eval(Expr(:using,:Gurobi))

    function EllTest1(cuts,flip)
        println(" Test1 cuts=$cuts, flip=$flip")
        m = RobustModel(solver=GurobiSolver(OutputFlag=0))
        @defVar(m, 0 <= x <= 10)
        @defUnc(m, 0 <= u <= 10)
        @setObjective(m, Max, 10x)
        !flip && addConstraint(m,  u*x <=  7)
         flip && addConstraint(m, -u*x >= -7)
        addEllipseConstraint(m, [1.0*u - 5], 2)  # 5 <= u <= 7
        #printRobust(m)
        solveRobust(m, suppress_warnings=true, prefer_cuts=cuts)
        #solveRobust(m,  prefer_cuts         = cuts,
        #                debug_printreform   =!cuts, 
        #                debug_printcut      = cuts,
        #                debug_printfinal    = true)
        @test_approx_eq_eps getValue(x) 1.0 1e-6
    end

    function EllTest2(cuts,flip)
        println(" Test2 cuts=$cuts, flip=$flip")
        m = RobustModel(solver=GurobiSolver(OutputFlag=0))
        @defVar(m, 0 <= x[i=1:5] <= 2*i)
        @defUnc(m, 0 <= u[i=1:5] <= i+4)
        @setObjective(m, Max, sum{(6-i)*x[i], i=1:5})
        !flip && addConstraint(m,  dot(u,x) <=  100)
         flip && addConstraint(m, -dot(u,x) >= -100)
        addEllipseConstraint(m, [3.0*u[1]-5, 1.0*u[5]-5, 2.0*u[4]-5], 1)
        #printRobust(m)
        solveRobust(m, suppress_warnings=true, prefer_cuts=cuts)
        @test_approx_eq_eps getValue(x[1]) 2.000 1e-6
        @test_approx_eq_eps getValue(x[2]) 4.000 1e-6
        @test_approx_eq_eps getValue(x[3]) 6.000 1e-6
        @test_approx_eq_eps getValue(x[4]) 8.000 1e-6
        @test_approx_eq_eps getValue(x[5]) 1.283 1e-2
    end

    function EllTest3(cuts,flip)
        # Polyhedral + ellipse
        println(" Test3 cuts=$cuts, flip=$flip")
        m = RobustModel(solver=GurobiSolver(OutputFlag=0))
        @defVar(m, 0 <= x[1:2] <= 10)
        @defUnc(m, u[1:2])
        @defUnc(m, z[1:2])
        @setObjective(m, Min, 1x[1] + 2x[2])
        !flip && addConstraint(m,  u[1]*x[1] + u[2]*x[2] >=  5)
         flip && addConstraint(m, -u[1]*x[1] - u[2]*x[2] <= -5)
        # Uncertainty set
        addConstraint(m, u[1] == 5.0*z[1]            + 10.0)
        addConstraint(m, u[2] == 3.0*z[1] - 2.0*z[2] +  3.0)
        addEllipseConstraint(m, [z[1],z[2]], 1)
        #printRobust(m)
        solveRobust(m, suppress_warnings=true, prefer_cuts=cuts)
        @test_approx_eq_eps getValue(x[1]) 1.000 1e-3
        @test_approx_eq_eps getValue(x[2]) 0.000 1e-3
    end

    function EllTest4(cuts,flip)
        println(" Test4 cuts=$cuts, flip=$flip")
        m = RobustModel(solver=GurobiSolver(OutputFlag=0))
        @defVar(m, 0 <= x <= 10)
        @defUnc(m, 0 <= u <= 10)
        @setObjective(m, Min, 10x)
        !flip && addConstraint(m,  x >=  u)
         flip && addConstraint(m, -x <= -u)
        addEllipseConstraint(m, [1.0*u - 5], 2)  # 5 <= u <= 7
        #printRobust(m)
        solveRobust(m, suppress_warnings=true, prefer_cuts=cuts) #, debug_printreform=!cuts, debug_printfinal=!cuts)
        @test_approx_eq_eps getValue(x) 7.0 1e-6
    end

    function EllTest5(cuts,flip)
        println(" Test5 cuts=$cuts, flip=$flip")
        m = RobustModel(solver=GurobiSolver(OutputFlag=0))
        @defVar(m, 0 <= x <=  8)
        @defUnc(m, 0 <= u <= 10)

        @defVar(m, 2 <= y <= 10)
        @defUnc(m, 0 <= w <= 10)

        @setObjective(m, Max, 20x + 10y)
        !flip && addConstraint(m,  u*x + w*y <=  10)
         flip && addConstraint(m, -u*x - w*y >= -10)
        addEllipseConstraint(m, [u - 5], 2)  # 5 <= u <= 7
        addEllipseConstraint(m, [w - 3], 1)  # 2 <= w <= 4
        #printRobust(m)
        solveRobust(m, suppress_warnings=true, prefer_cuts=cuts) #, debug_printreform=!cuts, debug_printfinal=!cuts)
        @test_approx_eq_eps getValue(x) (10-4*2)/7 1e-5
        @test_approx_eq_eps getValue(y) 2.0 1e-5
    end


    function EllTest6(cuts,flip)
        # Via @dfagnan
        println(" Test6 cuts=$cuts, flip=$flip")
        m = RobustModel(solver=GurobiSolver(OutputFlag=0))
        @defVar(m, 0 <= y <= 100)
        @defVar(m, 0 <= z <= 100)
        @defVar(m, -100 <= obj <= 100)
        @defUnc(m, 0 <= u[1:2] <= 100)
        @setObjective(m, Max, obj)
        if !flip
            addConstraint(m, obj <= (z+u[2]*y))
            addConstraint(m, z <= (1-u[1]*y))
        else
            addConstraint(m, -obj >= -(z+u[2]*y))
            addConstraint(m, -z >= -(1-u[1]*y))
        end
        addConstraint(m, u[1] == 1)
        addEllipseConstraint(m,[u[2]-1.2],0.01)
        #printRobust(m)
        solveRobust(m, suppress_warnings=true, prefer_cuts=cuts)
        @test_approx_eq_eps getValue(obj) 1.19 1e-6
    end

    for cuts in [true,false]
        for flip in [true,false]
            EllTest1(cuts,flip)
            EllTest2(cuts,flip)
            EllTest3(cuts,flip)
            EllTest4(cuts,flip)
            EllTest5(cuts,flip)
            EllTest6(cuts,flip)
        end
    end
end