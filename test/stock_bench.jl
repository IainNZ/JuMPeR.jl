using JuMP, JuMPeR
using Gurobi

r = 3.0
c = 1.0
h1 = 8.0
h2 = 4.0
M = 30.0
truck_cap = 5.0

immutable Sample
    D1
    D2
end

function solve_master(samples; silent=false, cuts=false)

    num_samp = length(samples)
    
    master = RobustModel(solver=GurobiSolver(OutputFlag=(silent?0:1), Threads=1))
    tic()

    @defVar(master, z)
    @defVar(master, zs[1:num_samp])
    @defVar(master, B[1:2] >= 0)
    @defVar(master, T[1:num_samp,1:2], Bin)
    @defVar(master, f[1:num_samp,1:2,0:2])  # S
    @defVar(master, g[1:num_samp,1:2,0:2])  # W

    @defUnc(master, D[1:num_samp,1:2])

    setObjective(master, :Max, 1.0*z)

    for samp_ind = 1:num_samp
        addConstraint(master, z <= zs[samp_ind])

        S1 = f[samp_ind,1,1]*D[samp_ind,1] + f[samp_ind,1,2]*D[samp_ind,2] + f[samp_ind,1,0]
        S2 = f[samp_ind,2,1]*D[samp_ind,1] + f[samp_ind,2,2]*D[samp_ind,2] + f[samp_ind,2,0]
        addConstraint(master, S1 >= 0)
        addConstraint(master, S2 >= 0)

        W1 = g[samp_ind,1,1]*D[samp_ind,1] + g[samp_ind,1,2]*D[samp_ind,2] + g[samp_ind,1,0]
        W2 = g[samp_ind,2,1]*D[samp_ind,1] + g[samp_ind,2,2]*D[samp_ind,2] + g[samp_ind,2,0]
        addConstraint(master, W1 >= 0)
        addConstraint(master, W2 >= 0)
        addConstraint(master, W1 <= truck_cap)
        addConstraint(master, W2 <= truck_cap)

        addConstraint(master, zs[samp_ind] <= r*(S1+S2) - c*(B[1]+B[2]) 
                                                - h1*T[samp_ind,1] - h2*T[samp_ind,2])
        addConstraint(master, S1 <= B[1] + W2 - W1)
        addConstraint(master, S2 <= B[2] + W1 - W2)
        addConstraint(master, S1 <= D[samp_ind,1])
        addConstraint(master, S2 <= D[samp_ind,2])
        addConstraint(master, W1 <= M*T[samp_ind,1])
        addConstraint(master, W2 <= M*T[samp_ind,2])

        addConstraint(master,  (D[samp_ind,1] - 20)/10 + (D[samp_ind,2] - 10)/5 <= 1)
        addConstraint(master,  (D[samp_ind,1] - 20)/10 - (D[samp_ind,2] - 10)/5 <= 1)
        addConstraint(master, -(D[samp_ind,1] - 20)/10 + (D[samp_ind,2] - 10)/5 <= 1)
        addConstraint(master, -(D[samp_ind,1] - 20)/10 - (D[samp_ind,2] - 10)/5 <= 1)

        for samp2 in 1:num_samp
            samp1 = samp_ind
            samp2 == samp1 && continue
            addConstraint(master, 
                                2*((samples[samp2].D1 - samples[samp1].D1)*D[samp_ind,1] +
                                   (samples[samp2].D2 - samples[samp1].D2)*D[samp_ind,2] )
                                    <=
                                samples[samp2].D1^2 + samples[samp2].D2^2 -
                                samples[samp1].D1^2 - samples[samp1].D2^2 )
        end
    end
    express_time = toq()
    !silent && println("Express time ", express_time)
    solveRobust(master, report=!silent, prefer_cuts=cuts)
end


samples = [Sample(20,10)]
srand(1998)
num_samples = 1
if length(ARGS) >= 1
    num_samples = int(ARGS[1])
end
println(num_samples)
for i = 1:(num_samples-1)
    D1 = rand(10:30)
    D2_max = +5*(1 - abs(D1-20)/10) + 10
    D2_min = -5*(1 - abs(D1-20)/10) + 10
    D2 = rand()*(D2_max-D2_min) + D2_min
    push!(samples, Sample(D1,D2))
end
println(samples)

cuts_flag = true

solve_master(samples, silent=true, cuts=cuts_flag)
using ProfileView
@profile solve_master(samples, silent=true, cuts=cuts_flag)
ProfileView.view()
readline()
for i = 1:5
    @time solve_master(samples, silent=true, cuts=cuts_flag)
end

