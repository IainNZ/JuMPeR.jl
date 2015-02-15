---
title: JuMPeR
author: Iain Dunning
order: 1
...

![JuMPeR Logo](logo.svg)\


# Introduction

**[JuMPeR]** is a modeling language for **[robust optimization (RO)]**. It is embedded in the **[Julia programming language]**, and is an extension to the **[JuMP]** modeling language. The goals of JuMPeR are  

* to make it easy to implement RO models,
* to make it easy to switch between solution methods, and
* to make it easy to try different uncertainty sets (via **oracles**.)

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
n, C = 10, 3
p, w = rand(n), rand(n)
σ = 0.2 * rand(n) .* w
knap = RobustModel()
@defVar(knap,      x[1:n], Bin)
@defUnc(knap, 0 <= z[1:n] <= 1)
@setObjective(knap, Max, sum{p[i]*x[i], i=1:n})
@addConstraint(knap, sum{(w[i]+σ[i]*z[i])*x[i], i=1:n} <= C)
solve(knap)
```

As JuMPeR builds on JuMP, you should be comfortable with the basics of modeling with [JuMP] before learning JuMPeR.

# Installation

JuMPeR is a registered Julia package, so you can simply install it with

```{.julia execute="false"}
Pkg.add("JuMPeR")
```

Note that both JuMP and JuMPeR do not come with a solver - you'll need to install one. You'll need to be a bit careful about making sure an appropriate solver is installed, for example:

* If we are solving a linear problem with an ellipsoidal uncertainty set, the reformulation will have second-order cone/quadratic constraints. For example, ECOS or Gurobi will work, but Clp will not.
* If we use a cutting-plane method to solve a problem with integer variables, the solver needs to support lazy constraints. For example GLPK and CPLEX will work, but Cbc will not.


# Modeling RO Problems

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




[Julia programming language]: http://julialang.org/
[JuMP]: https://github.com/JuliaOpt/JuMP.jl
[JuMPeR]: https://github.com/IainNZ/JuMPeR.jl
[robust optimization (RO)]: http://en.wikipedia.org/wiki/Robust_optimization