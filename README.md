JuMPeR
======
#### Julia for Mathematical Programming - extension for Robust optimization

[![Build Status](https://travis-ci.org/IainNZ/JuMPeR.jl.png?branch=master)](https://travis-ci.org/IainNZ/JuMPeR.jl)
[![Coverage Status](https://img.shields.io/coveralls/IainNZ/JuMPeR.jl.svg)](https://coveralls.io/r/IainNZ/JuMPeR.jl?branch=master)

**[JuMP]** is a domain-specific modeling language for **[mathematical programming]** embedded in **[Julia]**. It supports a number of open-source and commercial solvers.

**JuMPeR** extends JuMP to enable easy modeling of **[robust optimization]** problems. Specifically, instead of having to manually model the robust counterpart or write a cutting-plane generator, JuMPeR will do this work for you - simply provide the problem description, including uncertain parameters and the uncertainty set.

Documentation is available on [ReadTheDocs](http://jumper.readthedocs.org/en/latest/jumper.html)

[JuMP]: https://github.com/JuliaOpt/JuMP.jl
[mathematical programming]: http://en.wikipedia.org/wiki/Mathematical_optimization
[Julia]: http://julialang.org/
[robust optimization]: http://en.wikipedia.org/wiki/Robust_optimization

### Installation Instructions

JuMPeR isn't a listed package (yet). Heres what you're going to need to do to install it:

```julia
# You'll need JuMP, so install it if you haven't already
Pkg.add("JuMP")
# Now download JuMPeR direct from this repository
Pkg.clone("https://github.com/IainNZ/JuMPeR.jl.git")
# This will install it to your Julia package directory.
# Running Pkg.update() will always give you the freshest version of JuMPeR

# Of course, you'll need a solver, preferably one that supports lazy constraints.
# Gurobi is well tested with JuMPeR, but if you need a free one get GLPK:
Pkg.add("GLPKMathProgInterface")
```
