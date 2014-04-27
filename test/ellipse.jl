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