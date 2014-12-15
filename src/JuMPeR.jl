#############################################################################
# JuMPeR
# Julia for Mathematical Programming - extension for Robust Optimization
# See http://github.com/IainNZ/JuMPeR.jl
#############################################################################

module JuMPeR

# Bring all of JuMP into local namespace
importall JuMP  # Explicitly exported names
import JuMP.GenericAffExpr, JuMP.JuMPConstraint, JuMP.GenericRangeConstraint
import JuMP.sense, JuMP.rhs
import JuMP.IndexedVector, JuMP.addelt!, JuMP.isexpr
import JuMP.str_round
import JuMP.JuMPDict, JuMP.@gendict

export RobustModel, getNumUncs, solveRobust, printRobust
export setDefaultOracle!
export Uncertain, @defUnc, UAffExpr, FullAffExpr
export UncConstraint, UncSetConstraint

#############################################################################
# JuMP rexports
export
# Objects
    Model, Variable, AffExpr, QuadExpr, LinearConstraint, QuadConstraint,
    ConstraintRef, LinConstrRef,
# Functions
    # Model related
    getNumVars, getNumConstraints, getObjectiveValue, getObjective,
    getObjectiveSense, setObjectiveSense, writeLP, writeMPS, setObjective,
    addConstraint, addSOS1, addSOS2, solve,
    getInternalModel, buildInternalModel,
    # Variable
    setName, getName, setLower, setUpper, getLower, getUpper,
    getValue, setValue, getDual,
    # Expressions and constraints
    affToStr, quadToStr, conToStr, chgConstrRHS,
    
# Macros and support functions
    @addConstraint, @addConstraints, @defVar, 
    @defConstrRef, @setObjective, addToExpression, @defExpr, 
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
    uncCat::Vector{Symbol}

    # Adaptability
    adapt_type::Dict{Int,Symbol}
    adapt_on::Dict{Int,Vector}

    defaultOracle

    # Active cuts
    activecuts

    # Can have different solver for cutting planes
    cutsolver

    # For pretty printing
    dictList::Vector

    # Provided scenarios
    scenarios::Vector

    solve_time::Float64
end
RobustData(cutsolver) = RobustData(Any[],Any[],Any[],Any[],
                            0,String[],Float64[],Float64[],Symbol[],
                            Dict{Int,Symbol}(), Dict{Int,Vector}(),
                            GeneralOracle(), {},
                            cutsolver,{},{},0.0)

function RobustModel(;solver=JuMP.UnsetSolver(),cutsolver=JuMP.UnsetSolver())
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

getNumUncs(m::Model) = getRobust(m).numUncs

#############################################################################
# Uncertain
# Similar to JuMP.Variable, has an reference back to the model and an id num
type Uncertain
    m::Model
    unc::Int
end

function Uncertain(m::Model, lower::Number, upper::Number, cat::Symbol, name::String)
    robdata = getRobust(m)
    robdata.numUncs += 1
    push!(robdata.uncNames, name)
    push!(robdata.uncLower, convert(Float64, lower))
    push!(robdata.uncUpper, convert(Float64, upper))
    push!(robdata.uncCat, cat)
    return Uncertain(m, robdata.numUncs)
end
Uncertain(m::Model, lower::Number, upper::Number, cat::Symbol) = Uncertain(m,lower,upper,cat,"")

# Name setter/getters
setName(u::Uncertain, n::String) = (getRobust(u.m).uncNames[u.unc] = n)
function getName(u::Uncertain)
    checkUncNameStatus(u.m)
    return getRobust(u.m).uncNames[u.unc]
end
Base.print(io::IO, u::Uncertain) = print(io, getName(u))
Base.show( io::IO, u::Uncertain) = print(io, getName(u))

Base.isequal(u1::Uncertain, u2::Uncertain) = isequal(u1.unc, u2.unc)

#############################################################################
# Uncertain Affine Expression class
typealias UAffExpr GenericAffExpr{Float64,Uncertain}

UAffExpr() = UAffExpr(Uncertain[],Float64[],0.)
UAffExpr(c::Real) = UAffExpr(Uncertain[],Float64[],float(c))
UAffExpr(u::Uncertain) = UAffExpr([u],[1.],0.)
UAffExpr(u::Uncertain, c::Real) = UAffExpr([u],[float(c)],0.)
UAffExpr(coeffs::Array{Float64,1}) = [UAffExpr(c) for c in coeffs]
Base.zero(a::Type{UAffExpr}) = UAffExpr()  # For zeros(UAffExpr, dims...)
Base.zero(a::UAffExpr) = zero(typeof(a))

Base.print(io::IO, a::UAffExpr) = print(io, affToStr(a))
Base.show( io::IO, a::UAffExpr) = print(io, affToStr(a))


#############################################################################
# Full Affine Expression class
# Todo: better name. In my other robust modelling tools I called it
# something like this, but the catch then was that there we only two types of
# affexpr - the one with UAffExpr coefficients = Full, and the UAffExpr itself
typealias FullAffExpr GenericAffExpr{UAffExpr,Variable}

FullAffExpr() = FullAffExpr(Variable[], UAffExpr[], UAffExpr())
Base.zero(a::Type{FullAffExpr}) = FullAffExpr()
Base.zero(a::FullAffExpr) = zero(typeof(a))
function Base.push!(faff::FullAffExpr, new_coeff::Real, new_var::Variable)
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
    rd = getRobust(m)
    # Handle the odd special case where there are actually no variables in
    # the constraint - arises from use of macros
    if length(c.terms.vars) == 0
        # Pure uncertain
        @assert w == nothing
        return addConstraint(m, UncSetConstraint(c.terms.constant, c.lb, c.ub))
    end
    # Just a regular old constraint
    push!(rd.uncertainconstr,c)
    push!(rd.oracles, w)
    push!(rd.activecuts, {})
    return ConstraintRef{UncConstraint}(m,length(rd.uncertainconstr))
end

#############################################################################
# Scenarios
include("scenario.jl")

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

# Graph algorithms
include("graph.jl")

#############################################################################
end  # module
#############################################################################

#VG Talk to Iain about this
#Doesn't feel like best way
# Adding the DDUS stuff
include("../contrib/DDUS/ddusets.jl")

