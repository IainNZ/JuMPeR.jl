JuMPeR
======
#### Julia for Mathematical Programming - extension for Robust optimization

**JuMP** is a domain-specific modeling language for **[mathematical programming]**
embedded in **[Julia]**. It currently supports a number of open-source and
commercial solvers ([COIN Clp], [COIN Cbc], [GNU GLPK], [Gurobi], [MOSEK], [CPLEX]) via a 
[generic solver-independent interface](https://github.com/JuliaOpt/MathProgBase.jl). 

**JuMPeR** extends JuMP by making it possible to easily model robust optimization problems in a way that decouples the statement of the problem from both the specific uncertainty sets used and the solution method.

Build status: [![Build Status](https://travis-ci.org/IainNZ/JuMPeR.jl.png)](https://travis-ci.org/IainNZ/JuMPeR.jl)



[mathematical programming]: http://en.wikipedia.org/wiki/Mathematical_optimization
[Julia]: http://julialang.org/
[COIN Clp]: https://github.com/mlubin/Clp.jl
[COIN Cbc]: https://github.com/mlubin/Cbc.jl
[GNU GLPK]: http://www.gnu.org/software/glpk/
[Gurobi]: http://www.gurobi.com/
[MOSEK]: http://mosek.com/
[CPLEX]: http://www-01.ibm.com/software/commerce/optimization/cplex-optimizer/
