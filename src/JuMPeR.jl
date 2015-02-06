#############################################################################
# JuMPeR
# Julia for Mathematical Programming - extension for Robust Optimization
# See http://github.com/IainNZ/JuMPeR.jl
#############################################################################

module JuMPeR

# Import everything we need from JuMP, so we can build on it
importall JuMP
import JuMP.GenericAffExpr, JuMP.JuMPConstraint, JuMP.GenericRangeConstraint
import JuMP.sense, JuMP.rhs
import JuMP.IndexedVector, JuMP.addelt!, JuMP.isexpr
import JuMP.JuMPContainer, JuMP.JuMPDict, JuMP.JuMPArray
import JuMP.@gendict

export RobustModel, getNumUncs, solveRobust, printRobust
export setDefaultOracle!
export Uncertain, @defUnc, UAffExpr, FullAffExpr
export UncConstraint, UncSetConstraint



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
    setSolveHook
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
getName(u::Uncertain) = unc_str(REPLMode, u.m, u.unc)

Base.isequal(u1::Uncertain, u2::Uncertain) = isequal(u1.unc, u2.unc)

#############################################################################
# Uncertain Affine Expression class
typealias UAffExpr GenericAffExpr{Float64,Uncertain}

UAffExpr() = UAffExpr(Uncertain[],Float64[],0.)
UAffExpr(c::Real) = UAffExpr(Uncertain[],Float64[],float(c))
UAffExpr(u::Uncertain) = UAffExpr([u],[1.],0.)
UAffExpr(u::Uncertain, c::Real) = UAffExpr([u],[float(c)],0.)
# Next is a bit weird - its basically the vectorized version of UAffExpr(c) 
UAffExpr(coeffs::Array{Float64,1}) = [UAffExpr(c) for c in coeffs]
Base.zero(a::Type{UAffExpr}) = UAffExpr()  # For zeros(UAffExpr, dims...)
Base.zero(a::UAffExpr) = zero(typeof(a))

#############################################################################
# Full Affine Expression class
# Todo: better name. In my other robust modelling tools I called it
# something like this, but the catch then was that there we only two types of
# affexpr - the one with UAffExpr coefficients = Full, and the UAffExpr itself
typealias FullAffExpr GenericAffExpr{UAffExpr,Variable}

FullAffExpr() = FullAffExpr(Variable[], UAffExpr[], UAffExpr())
Base.zero(a::Type{FullAffExpr}) = FullAffExpr()
Base.zero(a::FullAffExpr) = zero(typeof(a))
function Base.push!(faff::FullAffExpr, new_coeff::UAffExpr, new_var::Variable)
    push!(faff.vars, new_var)
    push!(faff.coeffs, new_coeff)
end
function Base.push!(faff::FullAffExpr, new_coeff::Union(Real,Uncertain), new_var::Variable)
    push!(faff.vars, new_var)
    push!(faff.coeffs, UAffExpr(new_coeff))
end

#############################################################################
# UncSetConstraint      Just uncertainties
typealias UncSetConstraint GenericRangeConstraint{UAffExpr}
addConstraint(m::Model, c::UncSetConstraint) = push!(getRobust(m).uncertaintyset, c)

# UncConstraint         Mix of variables and uncertains
typealias UncConstraint GenericRangeConstraint{FullAffExpr}
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