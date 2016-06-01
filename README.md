![JuMPeR Logo](http://iainnz.github.io/JuMPeR.jl/logo.svg)

[![Build Status](https://travis-ci.org/IainNZ/JuMPeR.jl.svg?branch=master)](https://travis-ci.org/IainNZ/JuMPeR.jl)
[![codecov](https://codecov.io/gh/IainNZ/JuMPeR.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/IainNZ/JuMPeR.jl)

**[JuMPeR]** is a modeling language for [robust optimization].
It is embedded in the [Julia programming language], and is an extension to the [JuMP] modeling language.

* [Read documentation to learn more](doc/index.md)
* Install with Julia package manager: `julia> Pkg.add("JuMPeR")`

#### JuMPeR v0.4

JuMPeR was recently updated to version v0.4. This changed the syntax to match the changes in JuMP v0.13. Unlike JuMP v0.13, there are no deprecation warnings. However, it should generally be a simple find-and-replace, e.g. `@defUnc` to `@uncertain`. See the documentation for examples.

[Julia programming language]: http://julialang.org/
[JuMP]: https://github.com/JuliaOpt/JuMP.jl
[JuMPeR]: https://github.com/IainNZ/JuMPeR.jl
[robust optimization]: http://en.wikipedia.org/wiki/Robust_optimization
