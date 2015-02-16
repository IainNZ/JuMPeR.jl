---
title: JuMPeR
author: Iain Dunning
order: 1
...

![JuMPeR Logo](logo.svg)\


# About

**[JuMPeR]** is a modeling language for [robust optimization (RO)]. It is embedded in the [Julia programming language], and is an extension to the [JuMP] modeling language. JuMPeR was created by [Iain Dunning](http://iaindunning.com). Some of the goals of JuMPeR are

* to make it easy to implement RO models,
* to make it easy to switch between solution methods, and
* to make it easy to try different uncertainty sets (via 
<a href="#oracles">oracles</a>.)

Here is a small example of a robust 0/1 knapsack problem:
$$
\begin{alignat}{2}
\max \  & \sum_{i=1}^{n} p_i x_i \\
\text{subject to} \  & \sum_{i=1}^{n} \left( w_i + \hat{w}_i z_i \right) x_i \leq C \forall \mathbf{z} \in U \\
& \mathbf{x} \in \left\{0,1\right\}^n \\
\text{where} \  & U = \left\{ \mathbf{z} \in \left[0,1\right]^n, \sum_{i=1}^n z_i \leq \Gamma \right\}
\end{alignat}
$$
which we can write in JuMPeR as
```{.julia execute="false"}
using JuMP, JuMPeR
n, C, Γ = 10, 3, 4
p, w = rand(n), rand(n)
σ = 0.2 * rand(n) .* w
knap = RobustModel()
@defVar(knap,      x[1:n], Bin)
@defUnc(knap, 0 <= z[1:n] <= 1)
@setObjective(knap, Max, sum{p[i]*x[i], i=1:n})
@addConstraint(knap, sum{(w[i]+σ[i]*z[i])*x[i], i=1:n} <= C)
@addConstraint(knap, sum{z[i], i=1:n} <= Γ)
solve(knap)
```

As JuMPeR builds on JuMP, you should be comfortable with the basics of modeling with [JuMP] before learning JuMPeR.

# Installation

JuMPeR is a registered Julia package, so you can simply install it with

```{.julia execute="false"}
Pkg.add("JuMPeR")
```

Note that both JuMP and JuMPeR do not come with a solver - you'll need to install one. You'll need to be a bit careful about making sure an appropriate solver is installed, for example:

* If we are solving a linear problem with an ellipsoidal uncertainty set, the reformulation will have second-order cone/quadratic constraints. For example, [ECOS] or [Gurobi] will work, but [Clp] will not.
* If we use a cutting-plane method to solve a problem with integer variables, the solver needs to support lazy constraints. For example [GLPK] and [CPLEX] will work, but [Cbc] will not.

<a id="basics"></a>

# Basics of Modeling RO Problems

In this section we'll dive straight into modeling RO problems with the default *oracle*, and defer discussion of exactly what an oracle is (in the JuMPeR context) until later. You can think of them as a black box that JuMPeR asks how to solve a problem. They can reformulate constraints with uncertain parameters, and approximate them with cutting planes. The default oracle, called `GeneralOracle`, takes polyhedral and ellipsoidal constraints on the uncertain parameters and will reformulate them using duality or generate cutting planes by solving LPs/SOCPs.

## RobustModel

It all starts with the `RobustModel`, which is essentially the same as a JuMP `Model`. Like a `Model`, it also accepts a `solver` keyword argument. It can also take an additional `cutsolver` keyword argument, which the oracle can use for cutting planes.

```{.julia execute="false"}
# Use default solver for the problem class
rm = RobustModel()
# Can also specify the solver
using Gurobi
rm = RobustModel(solver=GurobiSolver())
# That'd produce a lot of output when using
# cutting planes, so we can set a solver just for cuts
rm = RobustModel(   solver=GurobiSolver(MIPGap=0.01),
                 cutsolver=GurobiSolver(OutputFlag=0))
```


## Uncertain Parameters

JuMPeR adds the `Uncertain` type to allow the user to model *uncertain parameters*. `Uncertain`s behave much like JuMP's `Variable`s - they have bounds and types, and can be combined with numbers, data, and inequalities to make expressions and constraints. They are defined with `@defUnc`:

```{.julia execute="false"}
rm = RobustModel()
# Can define single uncertain parameterss
@defUnc(rm, u)
# Can define indexed uncertain parameters too
@defUnc(rm, v[4:10,-2:2])
# They can have bounds too
@defUnc(rm, 0 <= w[i=1:10] <= i^2)
# You can also have types: integer (Int) and binary (Bin)
# Support for these will depend on the oracle. With the default
# oracle, they will only work with cutting plane mode, as
# reformulation is not possible in general for these types.
idxset = [:Apple, :Banana, :Orange]
@defUnc(rm, fruits[idxset], Bin)
```

There are some minor restrictions on `Uncertain` to be aware of:

- You cannot currently have products of uncertain parameters, e.g. `u*u` is not allowed.
- You cannot use uncertain parameters in the objective. If you have an uncertain objective, you should make an auxiliary variable to represent the objective and turn the objective into a contraint (e.g. an epigraph formulation):

$$
\begin{alignat}{2}
\max_{\mathbf{x} \in X} \min_{\mathbf{c} \in U} & \quad \mathbf{c} \cdot \mathbf{x}
\end{alignat}
$$

```{.julia execute="false"}
rm = RobustModel()
@defVar(rm, x[1:n])
@defUnc(rm, c[1:n])
# Create auxiliary for the objective
@defVar(rm, t)
# Epigraph form
@setObjective(rm, Max, t)
@addConstraint(rm, t <= sum{c[i]*x[i],i=1:n})
# ... constraints x in X ...
# ... constraints c in U ...
solve(rm)
```


## Expressions and Constraints

In JuMP we have only data and variables, with one type of expression and one type constraint (mostly). In JuMPeR we have three types of expression, each with are corresponding type of constraint:

1. Just data and variables (internally: `AffExpr`)
2. Just data and uncertain parameters (internally: `UAffExpr`)
3. Just data, uncertain parameters, and variables (internally: `FullAffExpr`)

Type 1 is just like in JuMP - these a deterministic constraints, and are JuMPeR doesn't do anything with them. JuMPeR treats Type 2 constraints as *uncertainty set constraints* - that is, they partially define the uncertainty set for this problem. They can be used by the oracles to perform reformulations or generate cutting planes. Type 3 constraints are uncertain constraints that must be feasible for all realizations of the uncertain parameters - these are what the oracles are dealing with to solve the problem. Here is an example of each, from the knapsack problem defined above:

```{.julia execute="false"}
w = rand(n)
σ = 0.2 * rand(n) .* w

knap = RobustModel()
@defVar(knap,      x[1:n], Bin)
@defUnc(knap, 0 <= z[1:n] <= 1)

# Type 1 constraint: only on the variables
#   "must take at least one item"
@addConstraint(knap, sum{x[i], i=1:n} >= 1)

# Type 2 constraint: only on the uncertain parameters
#   "only Γ parameters can be at their upper bound"
@addConstraint(knap, sum{z[i], i=1:n} <= Γ)

# Type 3 constraint: both variables and uncertain parameters
#   "the total uncertain weight of items is less than capacity"
@addConstraint(knap, sum{(w[i]+σ[i]*z[i])*x[i], i=1:n} <= C)
```

TODO: Ellipsodial constraints


## Solving

By default, JuMPeR will use duality to reformulate the uncertain problem into a deterministic robust counterpart. However you can signal that you want cuts with the `prefer_cuts` option. You can also get a summary report of how the model was solved and how long the different solve components took by setting the `report` option.

```{.julia execute="false"}
solve(rm)  # Will use default, which is reformulation
solve(rm, prefer_cuts=true)  # Use cuts if possible
solve(rm, report=true)  # Get some information about the solution process
```

Note that the solution methods available depends on the chosen oracle's capabilities. The default oracle supports both, but some uncertainty sets may only support one or the other.

<a id="oracles"></a>

# Oracles

JuMPeR's design is focussed around *oracles*. Oracles are responsible for taking RO problems, or parts of them, and transforming them into something solveable. That could be reformulating uncertain constraints into deterministic constraints, it could be generating cutting-planes, or something else entirely. Oracles are also intimately connected to uncertainty sets - for example, we can provide data to an oracle from which it can generate cutting-planes - an uncertainty set  never needs to explictly constructed. Finally oracles are interchangeable - you can obtain oracles from others or create your own to allow others to explore the performance of different sets and implementations.

In this section we will describe how to use oracles other than the default oracle, how to make an oracle, and some design considerations for oracles.

## Using oracles

JuMPeR currently comes with three oracles:

 * `GeneralOracle`, the default oracle. It takes an explicit polyhedral/ellipsoidal representation of the uncertainty set and can generate a reformulation or use cutting-planes.
 * `GeneralGraphOracle`, a variation on `GeneralOracle`. This oracle attempts to discover if the uncertain parameters actually belong to seperate, disjunct uncertainty sets. This allows it to generate smaller reformulations, and a seperator cut generator for each set if using cutting-planes.
 * `BertSimOracle`, implements the uncertainty set described in the 2004 paper *The Price of Robustness* by Bertsimas and Sim. Will generate cutting-planes efficiently using sorting instead of solving an LP. The `GeneralOracle` should be used if a reformulation is desired.

TODO: set default oracle, attaching oracles to constraints.

## Making an oracle

To make an oracle, we first need to understand what happens when you call `solve` on a `RobustModel`:

1. A new JuMP `Model`, referred to as the *master*, is created with the same variables as the original `RobustModel` and with all the deterministic constraints.
2. For each constraint, `registerConstraint` is called for the oracle associated with that constraint. Any constraints that don't have a oracle explicitly provided use the default.
3. `setup` is called for each oracle, giving them time to do any general setup shared across constraints. For example, it may take the dual of the uncertainty set in order to more efficiently reformulate multiple constraints.
4. Each oracle is then given a chance to reformulate its constraints (`generateReform`). It will return the number of constraints reformulated (which may be zero).
5. Start solving
 * If the problem is continuous, the master problem will be solved.
 * If the problem has integer variables, the MIP solver will be started.
6. Oracles get a chance to generate cutting-planes.
 * If the problem is continuous, the cutting-planes are added to the master problem and the master problem is resolved. If no new constraints are added, we terminate the solve.
 * If the problem has integer variables, the cutting-planes are added as lazy constraints only at integer solutions. The solve terminates when the MIP solver finds an optimal integer solution and no new constraints are added.

We can now specific how to make an oracle. An oracle should be a Julia type that is a subtype of `AbstractOracle`, i.e.

```{.julia execute="false"}
type MyNewOracle <: JuMPeR.AbstractOracle
    some_internal_state
end
```

The oracle must implement four methods, regardless of its functionality (defined in `src/oracle.jl`)

```{.julia execute="false"}
# registerConstraint
# Notifies the oracle that it is responsible for this constraint, and 
# passes any preferences provided via the solveRobust() command.
function registerConstraint(ab::AbstractOracle, rm::Model, ind::Int, prefs)

# setup
# Gives oracle time to do any setup it needs to do. Called after all
# constraints have been registered. Examples of work that could be done here
# include transforming the uncertainty set and generating a cutting plane
# model. Will NOT be called multiple times.
function setup(ab::AbstractOracle, rm::Model, prefs)

# generateReform
# Called before the main loop, adds anything it wants to the model. Returns
# number of constraints reformulated. If reformulation not supported or 
# desired, simply return 0
function generateReform(ab::AbstractOracle, master::Model, rm::Model, inds::Vector{Int})

# generateCut
# Called in the main loop every iteration/every time an integer solution is
# found. Returns a vector of constraints which are added to the problem by
# the main solve loop. Return an empty list if there are no constraints to add.
# The optional "active" argument will be called if the user wants to know
# the active scenarios at optimality. This is still experimental!
generateCut(ab::AbstractOracle, master::Model, rm::Model, inds::Vector{Int}, active=false)
```

## Oracle design advice

Yet to come.

[Julia programming language]: http://julialang.org/
[JuMP]: https://github.com/JuliaOpt/JuMP.jl
[JuMPeR]: https://github.com/IainNZ/JuMPeR.jl
[robust optimization (RO)]: http://en.wikipedia.org/wiki/Robust_optimization

[Cbc]: http://github.com/JuliaOpt/Cbc.jl
[Clp]: http://github.com/JuliaOpt/Clp.jl
[CPLEX]: http://github.com/JuliaOpt/CPLEX.jl
[ECOS]: http://github.com/JuliaOpt/ECOS.jl
[GLPK]: http://github.com/JuliaOpt/GLPK.jl
[Gurobi]: http://github.com/JuliaOpt/Gurobi.jl