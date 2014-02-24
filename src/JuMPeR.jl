#############################################################################
# JuMPeR
# Julia for Mathematical Programming - extension for Robust Optimization
# See http://github.com/IainNZ/JuMPeR.jl
#############################################################################

module JuMPeR

import JuMP.GenericAffExpr, JuMP.JuMPConstraint, JuMP.GenericRangeConstraint
import JuMP.sense, JuMP.rhs
import JuMP.IndexedVector, JuMP.addelt, JuMP.isexpr
importall JuMP  # What does this do exactly?

export RobustModel, Uncertain, UAffExpr, FullAffExpr, @defUnc, solveRobust

# JuMP rexports
export
# Objects
    Model, Variable, AffExpr, QuadExpr, LinearConstraint, QuadConstraint, MultivarDict,
# Functions
    # Relevant to all
    print,show,
    # Model related
    getNumVars, getNumConstraints, getObjectiveValue, getObjective,
    getObjectiveSense, setObjectiveSense, writeLP, writeMPS, setObjective,
    addConstraint, addVar, addVars, solve, copy,
    # Variable
    setName, getName, setLower, setUpper, getLower, getUpper,
    getValue, setValue, getDual,
    # Expressions and constraints
    affToStr, quadToStr, conToStr, chgConstrRHS,
    # Macros and support functions
    @addConstraint, @defVar, 
    @defConstrRef, @setObjective, addToExpression


#############################################################################
# RobustData contains all extensions to the base JuMP model type
type RobustData
    # Variable-Uncertain mixed constraints
    uncertainconstr
    # Oracles associated with each uncertainconstr
    oracles
    # Uncertain-only constraints
    uncertaintyset
    
    # Uncertainty data
    numUncs::Int
    uncNames::Vector{String}
    uncLower::Vector{Float64}
    uncUpper::Vector{Float64}
end
RobustData() = RobustData(Any[],Any[],Any[],0,String[],Float64[],Float64[])

function RobustModel(;solver=nothing)
    m = Model(solver=solver)
    m.ext[:Robust] = RobustData()
    return m
end

function getRobust(m::Model)
    if haskey(m.ext, :Robust)
        return m.ext[:Robust]
    else
        error("This functionality is only available for RobustModels")
    end
end


#############################################################################
# Uncertain
# Similar to JuMP.Variable, has an reference back to the model and an id num
type Uncertain
    m::Model
    unc::Int
end

function Uncertain(m::Model, lower::Number, upper::Number, name::String)
    robdata = getRobust(m)
    robdata.numUncs += 1
    push!(robdata.uncNames, name)
    push!(robdata.uncLower, convert(Float64, lower))
    push!(robdata.uncUpper, convert(Float64, upper))
    return Uncertain(m, robdata.numUncs)
end
Uncertain(m::Model, lower::Number, upper::Number) = Uncertain(m,lower,upper,"")

# Name setter/getters
setName(u::Uncertain, n::String) = (getRobust(u.m).uncNames[u.unc] = n)
function getName(u::Uncertain)
    n = getRobust(u.m).uncNames[u.unc]
    return n == "" ? string("_unc", u.unc) : n
end
print(io::IO, u::Uncertain) = print(io, getName(u))
show( io::IO, u::Uncertain) = print(io, getName(u))


#############################################################################
# Uncertain Affine Expression class
# Natural variant on GenericAffExpr, should probably use it but was created
# before... only difference is really uncs instead of vars

type UAffExpr
    uncs::Array{Uncertain,1}
    coeffs::Array{Float64,1}
    constant::Float64
end

UAffExpr() = UAffExpr(Uncertain[],Float64[],0.)
UAffExpr(c::Float64) = UAffExpr(Uncertain[],Float64[],c)
UAffExpr(u::Uncertain, c::Float64) = UAffExpr([u],[c],0.)
UAffExpr(coeffs::Array{Float64,1}) = [UAffExpr(c) for c in coeffs]
zero(::Type{UAffExpr}) = UAffExpr()  # For zeros(UAffExpr, dims...)

print(io::IO, a::UAffExpr) = print(io, affToStr(a))
show( io::IO, a::UAffExpr) = print(io, affToStr(a))


function affToStr(a::UAffExpr, showConstant=true)
    if length(a.uncs) == 0
        return string(a.constant)
    end

    # Get reference to robust part of model
    robdata = getRobust(a.uncs[1].m)

    # Collect like terms
    indvec = IndexedVector(Float64, robdata.numUncs)
    for ind in 1:length(a.uncs)
        addelt(indvec, a.uncs[ind].unc, a.coeffs[ind])
    end

    # Stringify the terms
    termStrings = Array(ASCIIString, length(a.uncs))
    numTerms = 0
    for i in 1:indvec.nnz
        idx = indvec.nzidx[i]
        numTerms += 1
        termStrings[numTerms] = "$(indvec.elts[idx]) $(robdata.uncNames[idx])"
    end

    # And then connect them up with +s
    ret = join(termStrings[1:numTerms], " + ")
    
    if abs(a.constant) >= 0.000001 && showConstant
        ret = string(ret," + ",a.constant)
    end
    return ret
end


#############################################################################
# Full Affine Expression class
# Todo: better name. In my other robust modelling tools I called it
# something like this, but the catch then was that there we only two types of
# affexpr - the one with UAffExpr coefficients = Full, and the UAffExpr itself
typealias FullAffExpr GenericAffExpr{UAffExpr,Variable}

FullAffExpr() = FullAffExpr(Variable[], UAffExpr[], UAffExpr())

# Pretty cool that this is almost the same as normal affExpr
function affToStr(a::FullAffExpr, showConstant=true)
    if length(a.vars) == 0
        return string(a.constant)
    end

    # Stringify the terms
    termStrings = Array(ASCIIString, length(a.vars))
    numTerms = 0
    for i in 1:length(a.vars)
        numTerms += 1
        termStrings[numTerms] = "($(affToStr(a.coeffs[i]))) $(getName(a.vars[i]))"
    end

    # And then connect them up with +s
    ret = join(termStrings[1:numTerms], " + ")
    
    # TODO(idunning): Think more carefully about this
    if showConstant
        ret = string(ret," + ",affToStr(a.constant))
    end
    return ret
end

#############################################################################
# UncSetConstraint class
# A constraint just involving uncertainties
typealias UncSetConstraint GenericRangeConstraint{UAffExpr}
addConstraint(m::Model, c::UncSetConstraint) = push!(getRobust(m).uncertaintyset, c)

#############################################################################
# UncConstraint class
# A mix of variables and uncertains
typealias UncConstraint GenericRangeConstraint{FullAffExpr}
function addConstraint(m::Model, c::UncConstraint, w=nothing)
    push!(getRobust(m).uncertainconstr,c)
    push!(getRobust(m).oracles, w)
end

#############################################################################
# Operator overloads
include("robustops.jl")

# All functions related to actual solution
include("solve.jl")

# Oracles... to be name changed
include("oracle.jl")

# Macros for more efficient generation
include("robustmacro.jl")

#############################################################################
end  # module
#############################################################################