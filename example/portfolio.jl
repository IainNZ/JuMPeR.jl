#############################################################################
# JuMPeR
# Julia for Mathematical Programming - extension for Robust Optimization
# See http://github.com/IainNZ/JuMPeR.jl
#############################################################################
# portfolio2.jl
# Builds a simple polyhedral uncertainty set using covariance information
# and uses the inbuilt PolyhedralOracle to solve it
#############################################################################

using JuMPeR
using Distributions
using Gurobi

const NUM_ASSET = 10

# Take in "box" and Gamma from the commandline. "Box" is how many standard
# deviations the "standard normals" in the uncertainty set can move from
# their means. Gamma is the total number of standard deviations away from
# mean you can have, ala Bertsimas Sim '04. See model definition below for
# more details.
if length(ARGS) != 2
    error("Expected two arguments, see code. Try julia portfolio.jl 1 2")
end
box   = int(ARGS[1])
Gamma = int(ARGS[2])

# Seed the RNG for reproducibilities sake
srand(10)


#############################################################################
# Simulate returns of the assets
# - num_samples is number of samples to take
# - Returns matrix, samples in rows, assets in columns
#############################################################################
function generate_data(num_samples)
    data = zeros(num_samples, NUM_ASSET)

    # Linking factors
    beta = [(i-1.)/NUM_ASSET for i = 1:NUM_ASSET] 

    for sample_ind = 1:num_samples
        # Common market factor, mean 3%, sd 5%, truncate at +- 3 sd
        z = rand(Normal(0.03, 0.05))
        z = max(z, 0.03 - 3*0.05)
        z = min(z, 0.03 + 3*0.05)

        for asset_ind = 1:NUM_ASSET
            # Idiosyncratic contribution, mean 0%, sd 5%, truncated at +- 3 sd
            asset = rand(Normal(0.00, 0.05))
            asset = max(asset, 0.00 - 3*0.05)
            asset = min(asset, 0.00 + 3*0.05)
            data[sample_ind, asset_ind] = beta[asset_ind] * z + asset
        end
    end

    return data
end


#############################################################################
# Solve the robust portfolio problem
#
# max   obj
# s.t.  sum(x_i      for i in 1:n) == 1
#       sum(r_i x_i  for i in 1:n) >= obj
#       x_i >= 0
# Uncertainty set:
#       r = A z + mean
#       y = |z| / box
#       sum(y_i for i in 1:n) <= Gamma
#       |z| <= box, 0 <= y <= 1
# where A is such that A*A' = Covariance matrix
#############################################################################
function solve_portfolio(past_returns, box, Gamma, pref_cuts, reporting)

    # Create covariance matrix and mean vector
    covar = cov(past_returns)
    means = mean(past_returns, 1)
    
    # Idea: multivariate normals can be described as
    # r = A * z + mu
    # where A*A^T = covariance matrix
    # So put a "standard" uncertainty set around z
    A = round(chol(covar),2)

    # Setup the robust optimization model
    m = RobustModel(solver=GurobiSolver(OutputFlag=0))

    # Variables
    @defVar(m, obj)  # Put objective as constraint
    @defVar(m, x[1:NUM_ASSET] >= 0)

    # Uncertainties
    @defUnc(m,         r[1:NUM_ASSET]       )  # The returns
    @defUnc(m, -box <= z[1:NUM_ASSET] <= box)  # The "standard normals"
    @defUnc(m,    0 <= y[1:NUM_ASSET] <= 1  )  # |z|/box

    @setObjective(m, Max, obj)

    # Portfolio constraint
    addConstraint(m, sum(x) == 1)

    # The objective constraint - uncertain
    addConstraint(m, dot(r, x) >= obj)

    # Build uncertainty set
    # First, link returns to the standard normals
    for asset_ind = 1:NUM_ASSET
        addConstraint(m, r[asset_ind] ==
            dot(A[asset_ind, :][:], z) + means[asset_ind])
    end
    # Then link absolute values to standard normals
    for asset_ind = 1:NUM_ASSET
        addConstraint(m, y[asset_ind] >= -z[asset_ind] / box)
        addConstraint(m, y[asset_ind] >=  z[asset_ind] / box)
    end
    # Finally, limit how much the standard normals can vary from means
    addConstraint(m, sum(y) <= Gamma)

    # Solve it, report statistics on number of cutting planes etc.
    reporting && println("Solving model...")

    # Uncomment the following lines to see model before printing
    #printRobust(m)
    #println("")

    solveRobust(m,  report=reporting, 
                    prefer_cuts=pref_cuts,
                    debug_printcut=false, #true,
                    debug_printfinal=false) #true)

    # Solution
    if reporting
        println("Solved, cuts prefered = $pref_cuts")
        println(getValue(x))
        println("Objective value: ", getValue(obj))
        println("")
    end

    return getValue(x), getValue(obj)
end


# Generate the simulated data
past_returns = generate_data(1000)
# Run once to warm start it
# solve_portfolio(past_returns, box, Gamma, true,  false)
# x, obj = solve_portfolio(past_returns, box, Gamma, false, false)
# Run again to see how fast it can go
# solve_portfolio(past_returns, box, Gamma, true,  true)
x, obj = solve_portfolio(past_returns, box, Gamma, true, true)

# Simulate some more data
NUM_FUTURE = 1000
future_returns = generate_data(NUM_FUTURE)
future_z = future_returns*x[:]
sort!(future_z)
println("Selected solution summary stats")
println("Minimum: ", future_z[1]*100)
println("10%:     ", future_z[int(NUM_FUTURE*0.1)]*100)
println("20%:     ", future_z[int(NUM_FUTURE*0.2)]*100)
println("30%:     ", future_z[int(NUM_FUTURE*0.3)]*100)
println("Mean:    ", mean(future_z)*100)
println("Maximum: ", future_z[end]*100)