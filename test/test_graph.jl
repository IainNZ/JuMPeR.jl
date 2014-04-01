using JuMPeR

m = RobustModel()
@defUnc(m, u[1:5,1:3])

for group = 1:5
    addConstraint(m, u[group,1] + group^2*u[group,2] == 1/group*u[group,3] )
end

printRobust(m)