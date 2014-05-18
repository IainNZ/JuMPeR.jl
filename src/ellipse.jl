#############################################################################
# JuMPeR
# Julia for Mathematical Programming - extension for Robust Optimization
# See http://github.com/IainNZ/JuMPeR.jl
#############################################################################
# ellipse.jl
# Ellipsoidal uncertainty set support
#############################################################################

export EllipseConstraint
export build_ellipse_constraint
export addEllipseConstraint

# EllipseConstraint
# Capture uncertainty set constraints of the form  || F u + g ||_2 <= Gamma
type EllipseConstraint <: JuMP.JuMPConstraint
    F::Array{Float64, 2}
    u::Array{Int, 1}
    g::Array{Float64, 1}
    Gamma::Float64
end

function printEll(rm::Model, ell::EllipseConstraint)
    rd = getRobust(rm)
    rows,cols = size(ell.F)
    row_strs = {}
    for r = 1:rows
        row_str = ""
        for c = 1:cols
            if ell.F[r,c] != 0.0
                row_str *= "$(ell.F[r,c]) $(rd.uncNames[ell.u[c]]) + "
            end
        end
        row_str *= "$(ell.g[r])"
        push!(row_strs, row_str)
    end
    max_len = maximum(map(length,row_strs))
    for r = 1:rows-1
        println("|| " * rpad(row_strs[r], max_len, " ") * " ||")
    end
    println("|| " * rpad(row_strs[rows], max_len, " ") * " || <= $(ell.Gamma)\n")
end

# build_ellipse_constraint
# Given || vec ||_2 <= Gamma, return an EllipseConstraint by expanding 
# `vec` out to its full `Fu+g` form.
function build_ellipse_constraint(vec::Vector, Gamma::Float64)
    num_terms = length(vec)

    # In the first pass we determine a unique set of uncertainties
    # present so we can allocate the correct size F and u
    unc_map  = Dict{Int,Int}()
    rev_map  = Dict{Int,Int}()
    num_uncs = 0
    for v in vec
        if typeof(v) <: UAffExpr
            for ind in 1:length(v.vars)
                unc = v.vars[ind].unc
                if !(unc in keys(unc_map))
                    num_uncs += 1
                    unc_map[unc] = num_uncs
                    rev_map[num_uncs] = unc
                end
            end
        elseif typeof(v) <: Uncertain
            if !(v.unc in keys(unc_map))
                num_uncs += 1
                unc_map[v.unc] = num_uncs
                rev_map[num_uncs] = v.unc
            end
        else
            error("Can only add norm constraints for simple functions of uncertainties.")
        end
    end
    u = [rev_map[i] for i=1:num_uncs]

    # Allocate memory
    F = zeros(num_terms, num_uncs)
    g = zeros(num_terms)
    for i in 1:num_terms
        if typeof(vec[i]) <: UAffExpr
            g[i] = vec[i].constant
            for ind in 1:length(vec[i].vars)
                unc = vec[i].vars[ind].unc
                F[i,unc_map[unc]] += vec[i].coeffs[ind]
            end
        elseif typeof(vec[i]) <: Uncertain
            F[i,unc_map[vec[i].unc]] += 1.0
        end
    end

    return EllipseConstraint(F,u,g,Gamma)
end

function addEllipseConstraint(m::Model, vec::Vector, Gamma::Real)
    push!(  getRobust(m).normconstraints,
            build_ellipse_constraint(vec, float(Gamma)) )
end