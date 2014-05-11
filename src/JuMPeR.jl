#############################################################################
# JuMPeR
# Julia for Mathematical Programming - extension for Robust Optimization
# See http://github.com/IainNZ/JuMPeR.jl
#############################################################################

module JuMPeR

import JuMP.GenericAffExpr, JuMP.JuMPConstraint, JuMP.GenericRangeConstraint
import JuMP.sense, JuMP.rhs
import JuMP.IndexedVector, JuMP.addelt!, JuMP.isexpr
import JuMP.string_intclamp
import JuMP.JuMPDict
importall JuMP

import Base.dot, Base.sum, Base.push!, Base.isequal

export RobustModel, Uncertain, UAffExpr, FullAffExpr, @defUnc, solveRobust
export UncConstraint, UncSetConstraint, printRobust
export setAdapt!
export setDefaultOracle!
export addScenario

#############################################################################
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
    addConstraint, addVar, addVars, addSOS1, addSOS2, solve, copy,
    getInternalModel,
    # Variable
    setName, getName, setLower, setUpper, getLower, getUpper,
    getValue, setValue, getDual,
    # Expressions and constraints
    affToStr, quadToStr, conToStr, chgConstrRHS,
    # Macros and support functions
    @addConstraint, @defVar, 
    @defConstrRef, @setObjective, addToExpression,
    @setNLObjective, @addNLConstraint


#############################################################################
# RobustData contains all extensions to the base JuMP model type
type RobustData
    # Variable-Uncertain mixed constraints
    uncertainconstr
    # Oracles associated with each uncertainconstr
    oracles
    # Uncertain-only constraints
    uncertaintyset
    normconstraints
    
    # Uncertainty data
    numUncs::Int
    uncNames::Vector{String}
    uncLower::Vector{Float64}
    uncUpper::Vector{Float64}
    uncCat::Vector{Int}

    # Adaptability
    adapt_type::Dict{Int,Symbol}
    adapt_on::Dict{Int,Vector}

    defaultOracle

    # Active cuts
    activecuts::Vector{Vector{Float64}}

    # Can have different solver for cutting planes
    cutsolver

    # For pretty printing
    dictList::Vector

    # Provided scenarios
    scenarios::Vector
end
RobustData(cutsolver) = RobustData(Any[],Any[],Any[],Any[],
                            0,String[],Float64[],Float64[],Int[],
                            Dict{Int,Symbol}(), Dict{Int,Vector}(),
                            GeneralOracle(), Vector{Float64}[],
                            cutsolver,JuMPDict[],Any[])

function RobustModel(;solver=nothing,cutsolver=nothing)
    m = Model(solver=solver)
    m.ext[:Robust] = RobustData(cutsolver)
    return m
end

function getRobust(m::Model)
    if haskey(m.ext, :Robust)
        return m.ext[:Robust]
    else
        error("This functionality is only available for RobustModels")
    end
end

function addScenario(m::Model, scen::Dict)
    robdata = getRobust(m)
    push!(robdata.scenarios, scen)
end


#############################################################################
# Uncertain
# Similar to JuMP.Variable, has an reference back to the model and an id num
type Uncertain
    m::Model
    unc::Int
end

function Uncertain(m::Model, lower::Number, upper::Number, cat::Int, name::String)
    robdata = getRobust(m)
    robdata.numUncs += 1
    push!(robdata.uncNames, name)
    push!(robdata.uncLower, convert(Float64, lower))
    push!(robdata.uncUpper, convert(Float64, upper))
    push!(robdata.uncCat, cat)
    return Uncertain(m, robdata.numUncs)
end
Uncertain(m::Model, lower::Number, upper::Number, cat::Int) = Uncertain(m,lower,upper,cat,"")

# Name setter/getters
setName(u::Uncertain, n::String) = (getRobust(u.m).uncNames[u.unc] = n)
function getName(u::Uncertain)
    checkUncNameStatus(u.m)
    return getRobust(u.m).uncNames[u.unc]
end
print(io::IO, u::Uncertain) = print(io, getName(u))
show( io::IO, u::Uncertain) = print(io, getName(u))

isequal(u1::Uncertain, u2::Uncertain) = (u1.unc == u2.unc)

#############################################################################
# Uncertain Affine Expression class
typealias UAffExpr GenericAffExpr{Float64,Uncertain}

UAffExpr() = UAffExpr(Uncertain[],Float64[],0.)
UAffExpr(c::Real) = UAffExpr(Uncertain[],Float64[],float(c))
UAffExpr(u::Uncertain) = UAffExpr([u],[1.],0.)
UAffExpr(u::Uncertain, c::Real) = UAffExpr([u],[float(c)],0.)
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
function push!(faff::FullAffExpr, new_coeff::Real, new_var::Variable)
    push!(faff.vars, new_var)
    push!(faff.coeffs, UAffExpr(new_coeff))
end

#############################################################################
# UncSetConstraint      Just uncertainties
typealias UncSetConstraint GenericRangeConstraint{UAffExpr}
# For 0.2
UncSetConstraint(uaff::UAffExpr,x::Float64,y::Int) = UncSetConstraint(uaff,x,float(y))
UncSetConstraint(uaff::UAffExpr,x::Int,y::Float64) = UncSetConstraint(uaff,float(x),x)
UncSetConstraint(uaff::UAffExpr,x::Int,y::Int) = UncSetConstraint(uaff,float(x),float(y))
addConstraint(m::Model, c::UncSetConstraint) = push!(getRobust(m).uncertaintyset, c)

# UncConstraint         Mix of variables and uncertains
typealias UncConstraint GenericRangeConstraint{FullAffExpr}
# For 0.2
UncConstraint(faff::FullAffExpr,x::Float64,y::Int) = UncConstraint(faff,x,float(y))
UncConstraint(faff::FullAffExpr,x::Int,y::Float64) = UncConstraint(faff,float(x),x)
UncConstraint(faff::FullAffExpr,x::Int,y::Int)     = UncConstraint(faff,float(x),float(y))
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
# Ellipsoidal uncertainty set support
include("ellipse.jl")

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