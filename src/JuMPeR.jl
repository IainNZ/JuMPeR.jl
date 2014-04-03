#############################################################################
# JuMPeR
# Julia for Mathematical Programming - extension for Robust Optimization
# See http://github.com/IainNZ/JuMPeR.jl
#############################################################################

module JuMPeR

import JuMP.GenericAffExpr, JuMP.JuMPConstraint, JuMP.GenericRangeConstraint
import JuMP.sense, JuMP.rhs
import JuMP.IndexedVector, JuMP.addelt, JuMP.isexpr
import JuMP.string_intclamp
importall JuMP

import Base.dot

export RobustModel, Uncertain, UAffExpr, FullAffExpr, @defUnc, solveRobust
export UncConstraint, UncSetConstraint, printRobust
export setAdapt!

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

    # Adaptability
    adapt_type::Dict{Int,Symbol}
    adapt_on::Dict{Int,Vector}

    defaultOracle

    # Active cuts
    activecuts::Vector{Vector{Float64}}
end
RobustData() = RobustData(  Any[],Any[],Any[],
                            0,String[],Float64[],Float64[],
                            Dict{Int,Symbol}(), Dict{Int,Vector}(),
                            PolyhedralOracle(), Vector{Float64}[Float64[]])

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
typealias UAffExpr GenericAffExpr{Float64,Uncertain}

UAffExpr() = UAffExpr(Uncertain[],Float64[],0.)
UAffExpr(c::Float64) = UAffExpr(Uncertain[],Float64[],c)
UAffExpr(u::Uncertain, c::Float64) = UAffExpr([u],[c],0.)
UAffExpr(coeffs::Array{Float64,1}) = [UAffExpr(c) for c in coeffs]
zero(::Type{UAffExpr}) = UAffExpr()  # For zeros(UAffExpr, dims...)

print(io::IO, a::UAffExpr) = print(io, affToStr(a))
show( io::IO, a::UAffExpr) = print(io, affToStr(a))


#############################################################################
# Full Affine Expression class
# Todo: better name. In my other robust modelling tools I called it
# something like this, but the catch then was that there we only two types of
# affexpr - the one with UAffExpr coefficients = Full, and the UAffExpr itself
typealias FullAffExpr GenericAffExpr{UAffExpr,Variable}

FullAffExpr() = FullAffExpr(Variable[], UAffExpr[], UAffExpr())

#############################################################################
# UncSetConstraint      Just uncertainties
typealias UncSetConstraint GenericRangeConstraint{UAffExpr}
addConstraint(m::Model, c::UncSetConstraint) = push!(getRobust(m).uncertaintyset, c)

# UncConstraint         Mix of variables and uncertains
typealias UncConstraint GenericRangeConstraint{FullAffExpr}
function addConstraint(m::Model, c::UncConstraint, w=nothing)
    push!(getRobust(m).uncertainconstr,c)
    push!(getRobust(m).oracles, w)
end

#############################################################################
# Adaptability
function setAdapt!(x::Variable, atype::Symbol, uncs::Vector)
    !(atype in [:Affine]) && error("Unrecognized adaptability type '$atype'")
    all_uncs = Uncertain[]
    add_to_list(u::Uncertain)           = (push!(all_uncs, u))
    add_to_list(u::JuMP.JuMPDict{Uncertain}) = (all_uncs=vcat(all_uncs,u.innerArray))
    add_to_list(u::Array{Uncertain})    = (all_uncs=vcat(all_uncs,vec(u)))
    add_to_list{T}(u::T)                = error("Can only depend on Uncertains (tried to adapt on $T)")
    map(add_to_list, uncs)
    rd = getRobust(x.m)
    rd.adapt_type[x.col] = atype
    rd.adapt_on[x.col]   = all_uncs
end
setAdapt!(x::JuMP.JuMPDict{Variable}, atype::Symbol, uncs::Vector) =
    map((v)->setAdapt!(v, atype, uncs), x.innerArray)
setAdapt!(x::Array{Variable}, atype::Symbol, uncs::Vector) =
    map((v)->setAdapt!(v, atype, uncs), x)

#############################################################################
# Operator overloads
include("robustops.jl")

# All functions related to actual solution
include("solve.jl")

# Oracles... to be name changed
include("oracle.jl")

# Macros for more efficient generation
include("robustmacro.jl")

# Pretty printing
include("print.jl")

#############################################################################
end  # module
#############################################################################