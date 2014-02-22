JuMPeR
======
#### Julia for Mathematical Programming - extension for Robust optimization

**JuMP** is a domain-specific modeling language for **[mathematical programming]**
embedded in **[Julia]**. It currently supports a number of open-source and
commercial solvers ([COIN Clp], [COIN Cbc], [GNU GLPK], [Gurobi], [MOSEK], [CPLEX]) via a 
[generic solver-independent interface](https://github.com/JuliaOpt/MathProgBase.jl). 

**JuMPeR** extends JuMP by making it possible to easily model robust optimization problems in a way that decouples the statement of the problem from both the specific uncertainty sets used and the solution method.