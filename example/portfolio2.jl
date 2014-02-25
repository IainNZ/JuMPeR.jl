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

const NUM_ASSET = 100

#############################################################################
# Simulate the data
#############################################################################
function generate_data(N)
    print("Generating data...")
    srand(10)

    beta = [(i-1.)/NUM_ASSET for i = 1:NUM_ASSET] # Linking factors
    data = zeros(N, NUM_ASSET)

    for sample_ind = 1:N
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

    println(" done")
    return data
end

#############################################################################
# Solve the robust portfolio problem
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
    #A = eye(NUM_ASSET)

    # Setup the robust optimization model
    m = RobustModel(solver=GurobiSolver(OutputFlag=0))

    # Variables
    @defVar(m, obj)  # Put objective as constraint
    @defVar(m, x[1:NUM_ASSET] >= 0)

    # Uncertainties
    @defUnc(m,         r[1:NUM_ASSET]       )
    @defUnc(m, -box <= z[1:NUM_ASSET] <= box)  # The "standard normals"
    @defUnc(m,    0 <= y[1:NUM_ASSET] <= 1  )  # |z|/box

    @setObjective(m, Max, obj)

    # Portfolio constraint
    addConstraint(m, sum([      x[i] for i=1:NUM_ASSET ])     == 1)

    # The objective constraint - uncertain
    addConstraint(m, sum([ r[i]*x[i] for i=1:NUM_ASSET ]) - obj >= 0)

    # Build uncertainty set
    # First, link returns to the standard normals
    for asset_ind = 1:NUM_ASSET
        addConstraint(m, r[asset_ind] == 
            sum([ A[asset_ind, j] * z[j] for j=1:NUM_ASSET]) + means[asset_ind] )
    end
    # Then link absolute values to standard normals
    for asset_ind = 1:NUM_ASSET
        addConstraint(m, y[asset_ind] >= -z[asset_ind] / box)
        addConstraint(m, y[asset_ind] >=  z[asset_ind] / box)
    end
    # Finally, limit how much the standard normals can vary from means
    addConstraint(m, sum([ y[j] for j=1:NUM_ASSET ]) <= Gamma)

    # Solve it, report statistics on number of cutting planes etc.
    reporting && println("Solving model...")

    # Uncomment the following lines to see model before printing
    # HACK until import better printing from JuMP/src/print.jl
    # for i in 1:NUM_ASSET
        # setName(r[i], "r[$i]")
        # setName(z[i], "z[$i]")
        # setName(y[i], "y[$i]")
    # end
    # printRobust(m)
    # println("")

    solveRobust(m,  report=reporting, 
                    prefer_cuts=pref_cuts,
                    debug_printcut=false, #true,
                    debug_printfinal=false) #true)

    # Solution
    if reporting
        println("Solved, cuts prefered = $pref_cuts")
        println(getValue(x))
        println(getValue(obj))
        println("")
    end
end

# Generate the simulated data
past_returns = generate_data(120)
# Run once to warm start it
solve_portfolio(past_returns, 1, 1, true,  false)
solve_portfolio(past_returns, 1, 1, false, false)
# Run again to see how fast it can go
solve_portfolio(past_returns, 1, 1, true,  true)
solve_portfolio(past_returns, 1, 1, false, true)