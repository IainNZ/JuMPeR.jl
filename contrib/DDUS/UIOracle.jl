###
# UI Oracle 
###
# Only supports cutting planes
using Optim
import JuMPeR: registerConstraint, setup, generateCut, generateReform

export UIOracle
export suppFcn

type UIOracle <: AbstractOracle
    lbounds::Vector{Float64}
    ubounds::Vector{Float64}
    data_sort::Matrix{Float64}
    log_eps::Float64

    # Cutting plane algorithm
    qL::Vector{Float64}
    qR::Vector{Float64}
    cut_tol::Float64  ##defaults to 1e-6

    # Other options
    debug_printcut::Bool
end

#preferred constructor
function UIOracle(data, lbounds, ubounds, eps, delta; cut_tol = 1e-6) 
    N, d = size(data)
    data_sort = sort_data_cols(data)
    Gamma = KSGamma(delta, N)
    qL, qR = gen_ql_qr(N, Gamma)
    UIOracle(vec(lbounds), vec(ubounds), data_sort, log(1/eps), qL, qR, 1e-6, false)
end

#returns bool, and ustar for degen case
is_degen(d, Gamma, log_eps) = d * log(1/Gamma) <= log_eps
degen_case(xs, lbounds::Vector{Float64}, ubounds::Vector{Float64}) = [xs[i] >= 0 ? ubounds[i] : lbounds[i] for i =1:length(xs)]

function gen_ql_qr(N::Int, Gamma)
    qL = zeros(Float64, N+2)
    qL[1] = Gamma
    qL[2:ifloor(N * (1-Gamma)) + 1] = 1/N
    qL[ifloor(N * (1-Gamma)) + 2] = 1-sum(qL)
    @assert (abs(sum(qL)-1) <= 1e-10) "QL not normalized $(sum(qL))"
    return qL, qL[N+2:-1:1]
end

function suppFcnUI(xs, data, lbounds, ubounds, log_eps, Gamma; cut_sense=:Max, lam_min=1e-8, lam_max=1e2, xtol=1e-8)
    data_sort = sort_data_cols(data)
    qL, qR = gen_ql_qr(size(data_sort, 1), Gamma)
    suppFcnUI(xs, data_sort, lbounds, ubounds, log_eps, qL, qR, cut_sense, lam_min, lam_max, xtol)
end

#returns zstar, ustar
function suppFcnUI(xs, data_sort, lbounds, ubounds, log_eps, qL, qR, cut_sense, lam_min=1e-8, lam_max = 1e2, xtol=1e-8)
    toggle = 1
    if cut_sense == :Min
        xs = copy(-xs)
        toggle = -1
    end

    const N = size(data_sort, 1)
    const d = size(data_sort, 2)
    const Gamma = qL[1]

    if is_degen(d, Gamma, log_eps)
        ustar = degen_case(xs, lbounds, ubounds)
        return toggle*dot(xs, ustar), ustar
    end

    #extend the data with the bounds
    const data = [ lbounds' ; data_sort ; ubounds' ]

    function f(lam::Float64)
        term = lam * log_eps
        term2 = 0.0
        for i =1:d
            #corrected for overflow
            if xs[i] > 0
                t = exp(xs[i] * (data[:, i] - ubounds[i]) / lam)
                t = dot(qR, t)
                term2 += lam * (xs[i] * ubounds[i]/lam + log(t))
            else
                t = exp(xs[i] * (data[:, i] - lbounds[i]) / lam)
                t = dot(qL, t)
                term2 += lam * (xs[i] * lbounds[i]/lam + log(t))
            end
            term += term2
        end
        term
    end
    res = Optim.optimize(f, lam_min, lam_max)
    !res.converged && error("Lambda linesearch did not converge")
    lamstar = res.minimum
    obj1 = res.f_minimum
    ustar = zeros(Float64, d)
    for i = 1:d
        if xs[i] >= 0
            qstar = qR .* exp(xs[i]*data[:, i]/lamstar)
        else
            qstar = qL .* exp(xs[i]*data[:, i]/lamstar)
        end
        qstar /= sum(qstar)
        ustar[i] = dot(qstar, data[:, i])
    end
    toggle*dot(ustar, xs), ustar
end

#preferred interface
suppFcn(xs, w::UIOracle, cut_sense) = suppFcnUI(xs, w.data_sort, w.lbounds, w.ubounds, w.log_eps, w.qL, w.qR, cut_sense)

# JuMPeR alerting us that we're handling this contraint
registerConstraint(w::UIOracle, rm::Model, ind::Int, prefs) = 
	! get(prefs, :prefer_cuts, true) && error("Only cutting plane supported")

function setup(w::UIOracle, rm::Model, prefs)
    # Extract preferences we care about
    w.debug_printcut = get(prefs, :debug_printcut, false)
    w.cut_tol        = get(prefs, :cut_tol, w.cut_tol)
    rd = JuMPeR.getRobust(rm)
    @assert (rd.numUncs == size(w.data_sort, 2)) "Num Uncertainties $(rd.numUncs) doesn't match columns in data $(size(w.data_sort, 2))"
    @assert (length(rd.uncertaintyset) == 0) "Auxiliary constraints on uncertainties not yet supported"
end

function generateCut(w::UIOracle, m::Model, rm::Model, inds::Vector{Int}, active=false)
    new_cons = {}
    rd = JuMPeR.getRobust(rm)
    for ind in inds
        con = JuMPeR.get_uncertain_constraint(rm, ind)
        cut_sense, xs, lhs_const = JuMPeR.build_cut_objective(rm, con, m.colVal)
        zstar, ustar = suppFcnUI(xs, w.data_sort, w.lbounds, w.ubounds, w.log_eps, w.qL, w.qR, cut_sense)
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
generateReform(w::UIOracle, m::Model, rm::Model, inds::Vector{Int}) = 0
