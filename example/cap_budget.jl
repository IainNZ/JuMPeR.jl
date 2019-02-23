#-----------------------------------------------------------------------
# JuMPeR  --  JuMP Extension for Robust Optimization
# http://github.com/IainNZ/JuMPeR.jl
#-----------------------------------------------------------------------
# Copyright (c) 2016: Iain Dunning
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
#-----------------------------------------------------------------------
# example/cap_budget.jl
# Apply the method from "Multistage Robust Mixed Integer Optimization
# with Adaptive Partitions" by Bertsimas and Dunning to the two-stage
# capital budgeting problem first described in Hanasusanto et. al.'s
# paper "Two-stage robust integer programming".
#
#   max  ∑ᵢ profitᵢ(ξ) [xᵢ + Θ yᵢ(ξ)]
#    st  ∑ᵢ costᵢ(ξ) [xᵢ + yᵢ(ξ)] ≤ B
#        xᵢ, yᵢ(ξ) ∈ {0,1}
#
# Requires a mixed-integer linear optimization problem solver.
#-----------------------------------------------------------------------

using JuMP, JuMPeR, LinearAlgebra

"""
    TreeScenario

Stores the values of "active" uncertain parameters, as well as the
associated tree structure described in the paper.
"""
mutable struct TreeScenario
    ξ::Vector{Float64}
    parent
    children::Vector
end
is_leaf(t::TreeScenario) = isempty(t.children)

"""
    generate_data(seed, N)

Generates a set of data for an instance of size `N` using the distributions
described in the paper by Hanasusanto, et. al., where `seed` is the random
seed to use to help with replication.
"""
function generate_data(seed, N)
    srand(seed)
    c0 = rand(N)*10  # Nominal costs
    r0 = c0 ./ 5     # Nominal profit
    B  = sum(c0)/2   # Budget
    θ = 0.8          # Profit discount factor for late decisions
    Φ = rand(N,4)    # Factor matrices ...
    Ψ = rand(N,4)    # ... for the uncertainty set
    for i in 1:N
        Φ[i,:] ./= sum(Φ[i,:])
        Ψ[i,:] ./= sum(Ψ[i,:])
    end
    return θ, B, Φ, Ψ, c0, r0
end

"""
    solve_partitioned_problem(N, θ, B, Φ, Ψ, c0, r0, scenarios)

Solve the two-stage problem of size `N` with one partition of the uncertainty
set for every leaf scenario in the `scenario_tree`. At optimality, grows the
tree given the new scenarios obtained.
"""
function solve_partitioned_problem(N::Int, θ::Float64, B::Float64,
                                   Φ::Matrix{Float64}, Ψ::Matrix{Float64},
                                   c0::Vector{Float64}, r0::Vector{Float64},
                                   scenario_tree::Vector{TreeScenario})
    # scenario_tree is a vector containing every scenario in the tree
    # We will have one partition for every leaf scenario in the tree
    leaf_scenarios = filter(is_leaf, scenario_tree)
    P = length(leaf_scenarios)  # Number of partitions
    # Initialize the RO model
    rm = RobustModel()
    # Decision variables: whether to pursue a project or not
    # First stage, here-and-now decision
    @variable(rm, x[1:N], Bin)
    # Second stage, wait-and-see decision
    # One set of variables per partition
    @variable(rm, y[1:P,1:N], Bin)
    # Define the uncertain parameters, which belong to a hypercube
    @uncertain(rm, -1 <= ξ[1:4] <= 1)
    # From the uncertain parameters we obtain the cost and profit vectors
    cost   = (1 + 0.5*Φ*ξ) .* c0
    profit = (1 + 0.5*Ψ*ξ) .* r0
    # The objective function will be the minimum of the objective function
    # across all the partitions. Put a default upper bound, just so we don't
    # start off unbounded if we are using a cutting plane method.
    @variable(rm, obj <= 10)
    @objective(rm, Max, obj)
    # For each partition...
    profit_con_refs = []
    budget_con_refs = []
    for p in 1:P
        # Define the partition uncertainty set. We define multiple hyperplanes
        # by walking up the scenario tree from this leaf.
        us = JuMPeR.BasicUncertaintySet()
        current_scenario = leaf_scenarios[p]
        parent_scenario = current_scenario.parent
        # We keep going until we hit the root of the tree, which is a scenario
        # that has no parent
        while parent_scenario != nothing
            for sibling_scenario in parent_scenario.children
                if current_scenario == sibling_scenario
                    continue  # Don't partition against ourself!
                end
                ξ_sub = sibling_scenario.ξ - current_scenario.ξ
                ξ_add = sibling_scenario.ξ + current_scenario.ξ
                @constraint(us, dot(ξ_sub, ξ) <= dot(ξ_sub,ξ_add)/2)
            end
            # Move up the scenario tree
            current_scenario  = parent_scenario
            parent_scenario = current_scenario.parent
        end
        # Constrain objective function for this partition
        profit_con_ref =
            @constraint(rm, obj <= sum{profit[i] * (x[i] + θ*y[p,i]), i in 1:N},
                                uncset=us)
        push!(profit_con_refs, profit_con_ref)
        # Constrain spending for this partition
        budget_con_ref =
            @constraint(rm, sum{cost[i] * (x[i] + y[p,i]), i in 1:N} <= B,
                                uncset=us)
        push!(budget_con_refs, budget_con_ref)
        # Can only choose each project once
        for i in 1:N
            @constraint(rm, x[i] + y[p,i] <= 1)
        end
    end
    # Solve, will use reformulation. We pass disable_cuts=true because
    # we want to be able to use solvers like CBC that don't support lazy
    # constraint callbacks.
    solve(rm, disable_cuts=true, active_scenarios=true)
    # Extend the scenario tree
    for p in 1:P
        # Extract the active uncertain parameter values
        profit_scen = getscenario(profit_con_refs[p])
        profit_scen_ξ = [uncvalue(profit_scen, ξ[i]) for i in 1:4]
        # Create a new child in the tree under this leaf
        profit_child = TreeScenario(profit_scen_ξ, leaf_scenarios[p], [])
        # Same for budget
        budget_scen = getscenario(budget_con_refs[p])
        budget_scen_ξ = [uncvalue(budget_scen, ξ[i]) for i in 1:4]
        budget_child = TreeScenario(budget_scen_ξ, leaf_scenarios[p], [])
        # Add to the tree
        push!(leaf_scenarios[p].children, profit_child)
        push!(leaf_scenarios[p].children, budget_child)
        push!(scenario_tree, profit_child)
        push!(scenario_tree, budget_child)
    end
    # Return the objective function value and the first-stage solution
    getobjectivevalue(rm), getvalue(x)
end

"""
    solve_problem(seed, N)

Solve the problem for the parameters corresponding to the given random seed
and instance size `N`. Prints results are each of five iterations.
"""
function solve_problem(seed, N)
    θ, B, Φ, Ψ, c0, r0 = generate_data(seed, N)
    # Start with no partitions (i.e., one scenario)
    scenario_tree = [ TreeScenario(zeros(4),nothing,[]) ]
    # Iteration 1
    z_iter1, x_iter1 = solve_partitioned_problem(N, θ, B, Φ, Ψ, c0, r0, scenario_tree)
    # Iteration 2
    z_iter2, x_iter2 = solve_partitioned_problem(N, θ, B, Φ, Ψ, c0, r0, scenario_tree)
    # Iteration 3
    z_iter3, x_iter3 = solve_partitioned_problem(N, θ, B, Φ, Ψ, c0, r0, scenario_tree)
    # Iteration 4
    z_iter4, x_iter4 = solve_partitioned_problem(N, θ, B, Φ, Ψ, c0, r0, scenario_tree)
    # Iteration 5
    z_iter5, x_iter5 = solve_partitioned_problem(N, θ, B, Φ, Ψ, c0, r0, scenario_tree)
    @show z_iter1
    println(round(Int, x_iter1))
    @show z_iter2
    println(round(Int, x_iter2))
    @show z_iter3
    println(round(Int, x_iter3))
    @show z_iter4
    println(round(Int, x_iter4))
    @show z_iter5
    println(round(Int, x_iter5))
    @show length(scenario_tree)
end

# By default, solve instance of size 5
solve_problem(106, 5)
