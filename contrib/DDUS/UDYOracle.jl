###
# UDY Oracle 
###
# At the moment only supports cutting planes
import JuMPeR: registerConstraint, setup, generateCut, generateReform
import JuMP.UnsetSolver

using JuMP, Mosek  #Right now SDPs require Mosek

export UDYOracle

type UDYOracle <: AbstractOracle
    eps_::Float64
    gamma1::Float64
    gamma2::Float64
    muhat::Vector{Float64}
    covhat::Array{Float64, 2}
    C::Array{Float64, 2}  #C'C = inv(covhat)

    # Cutting plane algorithm
    cut_tol::Float64  ##defaults to 1e-6
    cut_model::Model
    cut_vars::Vector{Variable}

    # Other options
    debug_printcut::Bool
end

#preferred interface
function UDYOracle(data, eps_, delta; numBoots=10000, cut_tol=1e-6, debug_printcut=false)
    muhat  = vec(mean(data, 1))
    covhat = cov(data)
    C = chol(inv(covhat))
    gamma1 = bootDY_mu(data, delta/2, numBoots)
    gamma2 = bootDY_sigma(data, delta/2, numBoots)
    UDYOracle(eps_, gamma1, gamma2, muhat, covhat, C, 
                cut_tol, Model(solver=MosekSolver()), Variable[], debug_printcut)
end

UDYOracle(muhat, covhat, gamma1, gamma2, eps_; cut_tol=1e-6, debug_printcut=false) = 
    UDYOracle(eps_, gamma1, gamma2, muhat, covhat, chol(inv(covhat)), 
                cut_tol, Model(solver=MosekSolver()), Variable[], debug_printcut)

# JuMPeR alerting us that we're handling this contraint
registerConstraint(w::UDYOracle, rm::Model, ind::Int, prefs) = 
    ! get(prefs, :prefer_cuts, true) && error("Only cutting plane supported")

#VG Make setup fcns consistent in terms of checks...
function setup(wDY::UDYOracle, rm::Model, prefs)
    wDY.cut_model = Model(solver=MosekSolver())

    # Extract preferences we care about
    wDY.debug_printcut = get(prefs, :debug_printcut, false)
    wDY.cut_tol        = get(prefs, :cut_tol, wDY.cut_tol)

    rd = JuMPeR.getRobust(rm)
    d = size(wDY.covhat, 1)
    @assert (rd.numUncs == d) "Num Uncertainties doesn't match columns in data"
    @assert (length(rd.uncertaintyset) == 0) #does not support additional cnsts on unctertainties for now

    #Assumes bounds were set directly on the robust contraints
    #build the cut solver
    @defMatrixVar(wDY.cut_model, u[d])
    @defMatrixVar(wDY.cut_model, w[d])
    @defMatrixVar(wDY.cut_model, m[d])
    @defVar(wDY.cut_model, 0 <= lam <=1/wDY.eps_)
    @defSDPVar(wDY.cut_model, A[d, d])
    @defSDPVar(wDY.cut_model, Ahat[d, d])

    addConstraint(wDY.cut_model, norm(wDY.C * w) <= sqrt(wDY.gamma1)*lam)

    for i = 1:d
        @addConstraint(wDY.cut_model, rd.uncLower[i] <= u[i] <= rd.uncUpper[i])
        @addConstraint(wDY.cut_model, rd.uncLower[i]*(lam-1) <= m[i] )
        @addConstraint(wDY.cut_model, m[i] <= rd.uncUpper[i]*(lam-1))

        @addConstraint(wDY.cut_model, 
            wDY.muhat[i]*lam  == u[i] + m[i] + w[i])
    end

    Z1 = [lam-1  m';
          m      A]
    addConstraint(wDY.cut_model, Z1 >= 0)

    Z2 = [1. u';
         u  Ahat ]
    addConstraint(wDY.cut_model, Z2 >= 0)

    #add an auxiliary variable for the big matrix expression
    @defMatrixVar(wDY.cut_model, Z3[d, d])
    sigbar = (wDY.gamma2 * wDY.covhat + wDY.muhat * wDY.muhat')
    for i = 1:d
        for j = 1:d
            @addConstraint(wDY.cut_model, Z3[i, j] == sigbar[i,j] * lam)
        end
    end

    temp = wDY.muhat * w'
    addConstraint(wDY.cut_model, Z3-A-Ahat-temp-temp' >=0)
    
    #VG Duplicate to make other things easier
    @defVar(wDY.cut_model, u2[1:d])
    for i = 1:d
        addConstraint(wDY.cut_model, u2[i] == u[i])
    end
    wDY.cut_vars = u2[:]
end

function generateCut(w::UDYOracle, m::Model, rm::Model, inds::Vector{Int}, active=false)
    rd = JuMPeR.getRobust(rm)
    master_sol = m.colVal
    new_cons = {}

    for ind in inds
        # Update the cutting plane problem's objective
        con = JuMPeR.get_uncertain_constraint(rm, ind)
        #VG This should move to sparse notation
        cut_sense, unc_obj_coeffs, lhs_const = 
                JuMPeR.build_cut_objective(rm, con, master_sol)
        
        setObjective(w.cut_model, cut_sense, sum([unc_obj_coeffs[ix] * w.cut_vars[ix] for ix=1:length(w.cut_vars)]))
        cut_solve_status = solve(w.cut_model, suppress_warnings=true)
        cut_solve_status != :Optimal && error("Cutting plane problem failed: $cut_solve_status")
        ustar = getValue(w.cut_vars)
        zstar = getObjectiveValue(w.cut_model)

        lhs_of_cut = zstar + lhs_const

        # SUBJECT TO CHANGE: active cut detection
        if active
            push!(rd.activecuts[ind], 
                JuMPeR.cut_to_scen(ustar, 
                    JuMPeR.check_cut_status(con, lhs_of_cut, w.cut_tol) == :Active))
            continue
        end

        # Check violation
        if JuMPeR.check_cut_status(con, lhs_of_cut, w.cut_tol) != :Violate
            w.debug_printcut && JuMPeR.debug_printcut(rm ,m,w,lhs_of_cut,con,nothing)
            continue  # No violation, no new cut
        end
        
        # Create and add the new constraint
        new_con = JuMPeR.build_certain_constraint(m, con, ustar)
        w.debug_printcut && JuMPeR.debug_printcut(rm, m, w, lhs_of_cut, con, new_con)
        push!(new_cons, new_con)
    end
    return new_cons
end

generateReform(w::UDYOracle, m::Model, rm::Model, inds::Vector{Int}) = 0
