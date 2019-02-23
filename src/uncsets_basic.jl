#-----------------------------------------------------------------------
# JuMPeR  --  JuMP Extension for Robust Optimization
# http://github.com/IainNZ/JuMPeR.jl
#-----------------------------------------------------------------------
# Copyright (c) 2016: Iain Dunning
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
#-----------------------------------------------------------------------
# src/uncsets_basic.jl
# The default uncertainty set, use explicitly provided constraints on
# the uncertain parameters to either build and repeatedly solve a
# cutting plane model, or use duality to generate reformulations of
# uncertain constraints.
# Included by src/uncsets.jl
#-----------------------------------------------------------------------


"""
    BasicUncertaintySet

The default uncertainty set, use explicitly provided constraints on the
uncertain parameters to either build and repeatedly solve a cutting plane
model, or use duality to generate reformulations of uncertain constraints.
"""
mutable struct BasicUncertaintySet <: AbstractUncertaintySet
    use_cuts::Bool

    # Explicitly provided constraints, other than the bounds
    linear_constraints::Vector{UncSetConstraint}
    norm_constraints::Vector{UncSetNormConstraint}

    # Cutting plane-specific options
    cut_model::Model
    cut_vars::Vector{Variable}
    cut_tol::Float64

    # Reformulation structure, generated in setup
    num_dualvar::Int
    dual_A::Vector{Vector{Tuple{Int,Float64}}}
    dual_objs::Vector{Float64}
    dual_vartype::Vector{Symbol}
    dual_contype::Vector{Symbol}
    dual_ell_rhs_idxs::Vector{Vector{Int}}
    dual_ell_lhs_idxs::Vector{Int}
    dual_l1_rhs_idxs::Vector{Vector{Int}}
    dual_l1_lhs_idxs::Vector{Int}
    dual_ω_idxs::Vector{Vector{Int}}
    dual_ω′_idxs::Vector{Vector{Int}}
end
BasicUncertaintySet() = BasicUncertaintySet(
                    false,  # use_cuts
                    UncSetConstraint[], UncSetNormConstraint[],
                    # Cutting plane
                    JuMP.Model(), JuMP.Variable[], 0.0,
                    # Reformulation
                    0, Vector{Tuple{Int,Float64}}[], Float64[], Symbol[], Symbol[],
                    Vector{Int}[], Int[],           # Reformulation, 2-norm
                    Vector{Int}[], Int[],           # Reformulation, 1-norm
                    Vector{Int}[], Vector{Int}[])   # Reformulation, ∞-norm


# Accept explicitly provided constraints
function JuMP.addconstraint(us::BasicUncertaintySet, c::UncSetConstraint)
    push!(us.linear_constraints, c)
    return ConstraintRef{BasicUncertaintySet,UncSetConstraint}(us, length(us.linear_constraints))
end
function JuMP.addconstraint(us::BasicUncertaintySet, c::UncSetNormConstraint)
    push!(us.norm_constraints, c)
    return ConstraintRef{BasicUncertaintySet,UncSetConstraint}(us, length(us.norm_constraints))
end


"""
    setup_set(BasicUncertaintySet, ...)

Generate the cutting plane model &| precompute the reformulation's structure.
"""
function setup_set(us::BasicUncertaintySet, rm::Model, idxs::Vector{Int},
                    scens_requested::Bool, other_prefs::Dict{Symbol,Any})
    # Extract preferences we care about
    us.use_cuts = get(other_prefs, :prefer_cuts, false)
    us.cut_tol  = get(other_prefs, :cut_tol, 1e-6)
    # Set up only what is needed
    if us.use_cuts || scens_requested
        setup_set_cut(us, rm)
    end
    if !us.use_cuts
        setup_set_reform(us, rm)
    end
end

# Cutting plane-specific functionality
include("uncsets_basic_cut.jl")

# Reformuation-specifc functionality
include("uncsets_basic_reform.jl")
