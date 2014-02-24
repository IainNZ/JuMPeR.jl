JuMPeR
======
#### Julia for Mathematical Programming - extension for Robust optimization

**JuMP** is a domain-specific modeling language for **[mathematical programming]**
embedded in **[Julia]**. It currently supports a number of open-source and
commercial solvers ([COIN Clp], [COIN Cbc], [GNU GLPK], [Gurobi], [MOSEK], [CPLEX]) via a 
[generic solver-independent interface](https://github.com/JuliaOpt/MathProgBase.jl). 

**JuMPeR** extends JuMP by making it possible to easily model robust optimization problems in a way that decouples the statement of the problem from both the specific uncertainty sets used and the solution method.

Build status: [![Build Status](https://travis-ci.org/IainNZ/JuMPeR.jl.png)](https://travis-ci.org/IainNZ/JuMPeR.jl)

### Installation Instructions

JuMPeR isn't a listed package (yet). Heres what you're going to need to do to install it:

```
Pkg.add("JuMP")             # JuMP is a dependency
Pkg.checkout("JuMP")        # For now we are using cutting edge JuMP, but
                            # will probably be fine using normal after JuMP 0.4
Pkg.clone("https://github.com/IainNZ/JuMPeR.jl.git")
                            # Installs JuMPeR straight from the repository -
                            # now you can get updates using Pkg.update()
# If you don't have a solver, get a free one.
Pkg.add("Cbc")
Pkg.add("Clp")
# Optional other packages
Pkg.add("Distributions")    # Used by the portfolio example to generate data.
```


[mathematical programming]: http://en.wikipedia.org/wiki/Mathematical_optimization
[Julia]: http://julialang.org/
[COIN Clp]: https://github.com/mlubin/Clp.jl
[COIN Cbc]: https://github.com/mlubin/Cbc.jl
[GNU GLPK]: http://www.gnu.org/software/glpk/
[Gurobi]: http://www.gurobi.com/
[MOSEK]: http://mosek.com/
[CPLEX]: http://www-01.ibm.com/software/commerce/optimization/cplex-optimizer/
