![JuMPeR Logo](logo.svg)


**[JuMPeR](https://github.com/IainNZ/JuMPeR.jl)** is a modeling language for robust optimization (RO) (see, e.g., ["The Price of Robustness"](http://www.mit.edu/~dbertsim/papers/Robust%20Optimization/The%20price%20of%20Robustness.pdf)). It is embedded in the [Julia programming language](http://julialang.org/), and is an extension to the [JuMP](https://github.com/JuliaOpt/JuMP.jl) modeling language. JuMPeR was created by [Iain Dunning](http://iaindunning.com). JuMPeR supports the following features:

* Easily create polyhedral and ellipsoidal uncertainty sets
* Use a _reformulation_ or _cutting plane_ approach to solving the problems
* Special support for _affine adaptability_ (linear decision rules)
* Extensible uncertainty set system for more exotic uncertainty sets.

JuMPeR does't have a corresponding paper at this time. I'd appreciate if, for now, that you cite the source code and the [JuMP paper](http://arxiv.org/abs/1508.01982). As JuMPeR builds on JuMP, you should be comfortable with the basics of modeling with JuMP before learning JuMPeR.


# Installation

JuMPeR is a registered Julia package, so you can simply install it with

```jl
Pkg.add("JuMPeR")
```

Note that both JuMP and JuMPeR do not come with a solver - you'll need to install one listed on the [JuliaOpt website](http://juliaopt.org). Make sure that you have the required solver for your problem class, or you might receive an error. For example:

* For a linear problem with an ellipsoidal uncertainty set, the reformulation will have second-order cone/quadratic constraints. Solvers like `ECOS` or `Gurobi` will work, but solvers like `Clp` will not.
* If using cutting planes for a problem with integer variables, the solver must support _lazy constraints_. For example `GLPK` and `Gurobi` will work, but `Cbc` will not.


# JuMPeR Basics

JuMPeR introduces a new primitive for _uncertain parameters_, `Uncertain`. To use `Uncertain`s in your model, you must use a `RobustModel` instead of JuMP's `Model`, e.g.,

```jl
using JuMPeR
rm = RobustModel()
# Set the solver just like you would for Model
rm = RobustModel(solver=GurobiSolver(OutputFlag=0))
```

Uncertain parameters are defined with `@uncertain`, which behaves just like `@variable` (except that providing starting values and some categories are not supported):

```jl
rm = RobustModel()
# A single uncertain parameter
@uncertain(rm, capacity)
# A set of uncertain parameters, with bounds
@uncertain(rm, μ[i]-σ[i] <= weights[i=1:N] <= μ[i]+σ[i])
# You can even define binary uncertain parameters - although you might not
# be able to solve the resulting problem with all possible some methods!
# Also, note the use of an arbitrary index set.
highways = [:I90, :I93, :I95]
@uncertain(rm, blocked[highways], Bin)
```

You can use uncertain parameters in expressions much like variables, but with a couple of restrictions. The most important is that uncertain parameters are not supported in objective functions. If you have an uncertain objective function, you should use an _epigraph formulation_ by turning the objective function into a constraint:

```jl
rm = RobustModel()
@variable(rm, x);  @uncertain(rm, u)
# BAD:
@objective(rm, Max, u*x)
# Good:
@variable(rm, z)
@objective(rm, Max, z)
@variable(rm, z <= u*x)
```

The other restrictions are due to limits of JuMPeR at this time. In particular, JuMPeR doesn't support products of uncertain parameters (e.g. `u²`), and doesn't support using uncertain parameters in quadratic or second-order cone constraints. Everything else is fair game! To put it together, here is an example of a knapsack problem with uncertain item weights (and a box uncertainty set):

```jl
knap = RobustModel()
@variable(knap, take[1:N], Bin)
@uncertain(knap, weight_low[i] <= weight[1:N] <= weight_hi[i])
@objective(knap, Max, sum{profit[i]*take[i], i=1:N})
@constraint(knap, sum{weight[i]*take[i], i=1:N} <= capacity)
solve(knap)
```


# Working with Uncertainty Sets

Mathematically speaking, an uncertainty set describes the values the uncertain parameters can take in a constraint. In JuMPeR, an `UncertaintySet` is a user-defined type that controls how a problem is solved. JuMPeR comes with a `BasicUncertaintySet` that is sufficient for most users needs. We'll describe how the `UncertaintySet` interface works in SectionX.

When a `RobustModel` is created, JuMPeR creates a problem-wide default `BasicUncertaintySet` that will be used for all constraints. The `BasicUncertaintySet` allows us to add linear and norm constraints to the uncertain parameters -- constraints that only include uncertain parameters are automatically recognized as uncertainty set-defining constraints.

Consider a robust portfolio problem. We have an uncertain return $r_i$ for each asset $i$, which we'll model two (and a half) different ways. We'll first consider an ellipsoidal uncertainty set

${\cal U}_E = \left\{ \left(r,\xi\right) \mid
r = \Sigma \xi + \mu,
\ \Vert\xi\Vert_2 \leq \Gamma
\right\},$

where $\xi$ is an auxiliary uncertain parameter that captures how far we deviate from the nominal returns $\mu$, and $\Sigma$ is a covariance matrix that links the assets together. We can write this in JuMPeR with a combination of linear and norm constraints:

```jl
port = RobustModel()
@uncertain(port, r[1:N])
@uncertain(port, ξ[1:N])
@constraint(port, r .== Σ*ξ + μ)  # element-wise equality
@constraint(port, norm(ξ, 2) <= Γ)
```

The second set we'll consider is the "budget" polyhedral uncertainty set from the ["Price of Robustness" paper by Bertsimas and Sim](http://www.mit.edu/~dbertsim/papers/Robust%20Optimization/The%20price%20of%20Robustness.pdf)). In words, we want to allow no more than $\Gamma$ of the uncertain returns to deviate from their nominal values. We can express this mathematically in a form similar to the ellipsoidal set above, using norms:

${\cal U}_P = \left\{ \left(r,\xi\right) \mid
r = \Sigma \xi + \mu,
\ \Vert\xi\Vert_1 \leq \Gamma,
\ \Vert\xi\Vert_\infty \leq 1
\right\},$

We can approach this in two ways in JuMPeR. The first is to use norms, like above:

```jl
port = RobustModel()
@uncertain(port, r[1:N])
@uncertain(port, ξ[1:N])
@constraint(port, r .== Σ*ξ + μ)  # element-wise equality
@constraint(port, norm(ξ, Inf) <= 1)
@constraint(port, norm(ξ, 1) <= Γ)
```
