#####
# UM
#####
import JuMPeR: registerConstraint, setup, generateCut, generateReform
import JuMP.UnsetSolver

export UMOracle
export suppFcn

#returns zstar, ustar
function suppFcnUM(xs, lquants, uquants, cut_sense=:Max)
    const d = length(lquants)
    if cut_sense == :Min
    	ustar = [xs[i] > 0 ? lquants[i] : uquants[i] for i =1:d]
        return sum(min(lquants .* xs, uquants .* xs)), ustar
    else
    	ustar = [xs[i] > 0 ? uquants[i] : lquants[i] for i =1:d]
    	return sum(max(lquants .* xs, uquants .* xs)), ustar
    end
end

function calc_s(data, eps_, delta)
	N, d = size(data)
	if (1-eps_/d)^N  > delta / 2d
		return N + 1
	else
		dBin = Binomial(N, 1-eps_/d)
		return quantile(dBin, 1-delta/2d)
	end
end

######################
type UMOracle <: AbstractOracle
    lquants::Vector{Float64}
    uquants::Vector{Float64}
    cut_tol::Float64  ##defaults to 1e-6

    # Other options
    debug_printcut::Bool
end

suppFcn(xs::Vector, w::UMOracle, cut_sense) = 
        suppFcnUM(xs, w.lquants, w.uquants, cut_sense)

#Preferred Interface
function UMOracle(data, lbounds, ubounds, eps_, delta; cut_tol=1e-6, debug_printcut=false)
	s = calc_s(data, eps_, delta)
	N = size(data, 1)
	@assert N- s + 1 < s "UM: N not sufficiently big N: $N \t s: $s"
	if s == N + 1
		return UMOracle(lbounds, ubounds, cut_tol, debug_printcut)
	else
		data_sort = sort_data_cols(data)
		return UMOracle(vec(data_sort[N-s+1, :]), vec(data_sort[s, :]), cut_tol, debug_printcut)
	end
end

# JuMPeR alerting us that we're handling this contraint
registerConstraint(w::UMOracle, rm::Model, ind::Int, prefs) = 
    ! get(prefs, :prefer_cuts, true) && error("Only cutting plane supported")

function setup(w::UMOracle, rm::Model, prefs)
    # Extract preferences we care about
    w.debug_printcut = get(prefs, :debug_printcut, false)
    w.cut_tol        = get(prefs, :cut_tol, w.cut_tol)
    rd = JuMPeR.getRobust(rm)
    @assert (rd.numUncs == length(w.lquants)) "Num Uncertainties $(rd.numUncs) doesn't match columns in lquants $(w.lquants)"
    @assert (rd.numUncs == length(w.uquants)) "Num Uncertainties $(rd.numUncs) doesn't match columns in uquants $(w.uquants)"
    @assert (length(rd.uncertaintyset) == 0) "Auxiliary constraints on uncertainties not yet supported"
end

#Refactor all this logic into a single place
function generateCut(w::UMOracle, m::Model, rm::Model, inds::Vector{Int}, active=false)
    new_cons = {}
    rd = JuMPeR.getRobust(rm)
    for ind in inds
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
generateReform(w::UMOracle, m::Model, rm::Model, inds::Vector{Int}) = 0






