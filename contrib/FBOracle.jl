###
# FB Oracle 
###
# At the moment only supports cutting planes
import JuMPeR: registerConstraint, setup, generateCut, generateReform

export FBOracle
export suppFcn

type FBOracle <: AbstractOracle
    mfs::Vector{Float64}
    mbs::Vector{Float64}
    sigfs::Vector{Float64}
    sigbs::Vector{Float64}
    log_eps::Float64

    # Cutting plane algorithm
    cut_tol::Float64  ##defaults to 1e-6

    # Other options
    debug_printcut::Bool
end

FBOracle(mfs, mbs, sigfs, sigbs, eps; TOL=1e-6) = 
    FBOracle(mfs, mbs, sigfs, sigbs, log(1/eps), TOL, false)

#Preferred constructors
function FBOracle(data, eps, delta1, delta2; CUT_TOL=1e-6, numBoots=int(1e4))
    N, d  = size(data)
    mfs   = zeros(Float64, d)
    mbs   = zeros(Float64, d)
    sigfs = zeros(Float64, d)
    sigbs = zeros(Float64, d)

    for i = 1:d
        mbs[i], mfs[i] = calcMeansT(data[:, i], delta1/d)
        sigfs[i], sigbs[i] = calcSigsBoot(data[:, i], delta2/d, numBoots)
    end
    FBOracle(mfs, mbs, sigfs, sigbs, eps, TOL=CUT_TOL)
end
FBOracle(data, eps, delta; CUT_TOL=1e-6, numBoots=int(1e4)) = 
    FBOracle(data, eps, delta/2, delta/2, CUT_TOL=CUT_TOL, numBoots=numBoots)

suppFcn(xs, w::FBOracle, cut_sense) = 
    suppFcnFB(xs, w.mfs, w.mbs, w.sigfs, w.sigbs, w.log_eps, cut_sense)

#log_eps = log(1/eps_)
#returns zstar, ustar
function suppFcnFB(xs, mfs, mbs, sigfs, sigbs, log_eps, cut_sense=:Max)
    sign_flip = 1
    if cut_sense == :Min
        xs = copy(-xs)
        sign_flip = -1
    end
    lam = 0
    for i = 1:length(xs)
        if xs[i] >= 0
            lam += sigfs[i]^2 * xs[i]^2
        else
            lam += sigbs[i]^2 * xs[i]^2
        end
    end
    lam /= (2 * log_eps)
    lam = sqrt(lam)
    ustar = zeros(length(xs))
    for i = 1:length(xs)
        if xs[i] >= 0
            ustar[i] = mfs[i] + xs[i] * sigfs[i]^2/lam
        else
            ustar[i] = mbs[i] + xs[i] * sigbs[i]^2/lam
        end
    end
    zstar = dot(xs, ustar) * sign_flip
    return zstar, ustar
end



# JuMPeR alerting us that we're handling this contraint
registerConstraint(w::FBOracle, rm::Model, ind::Int, prefs) = 
	! get(prefs, :prefer_cuts, true) && error("Only cutting plane supported")

function setup(w::FBOracle, rm::Model, prefs)
    # Extract preferences we care about
    w.debug_printcut = get(prefs, :debug_printcut, false)
    w.cut_tol        = get(prefs, :cut_tol, w.cut_tol)

    rd = JuMPeR.getRobust(rm)
    @assert (rd.numUncs == length(w.mfs)) "Num Uncertainties $(rd.numUncs) doesn't match columns in data $(size(w.mfs))"
    @assert (length(w.mfs) == length(w.mbs) == length(w.sigfs) == length(w.sigbs)) "Lengths of means and fb devs dont match uncertainties"

    #ignore any additional constraints on uncertainties for now
    @assert (length(rd.uncertaintyset) == 0) "Auxiliary constraints on uncertainties not yet supported"
end

function generateCut(w::FBOracle, m::Model, rm::Model, inds::Vector{Int}, active=false)
    new_cons = {}
    rd = JuMPeR.getRobust(rm)

    for ind in inds
        # Update the cutting plane problem's objective
        con = JuMPeR.get_uncertain_constraint(rm, ind)
        cut_sense, xs, lhs_const = JuMPeR.build_cut_objective(rm, con, m.colVal)
        zstar, ustar = suppFcn(xs, w, cut_sense) 
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

#Shouldn't be called
generateReform(w::FBOracle, m::Model, rm::Model, inds::Vector{Int}) = 0


