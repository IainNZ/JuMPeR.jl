###
# LCX Oracle 
###
# At the moment only supports cutting planes
using JuMPeR
import JuMPeR: registerConstraint, setup, generateCut, generateReform
import MathProgBase.AbstractMathProgSolver

export LCXOracle
export suppFcn

#returns bool, violated ab_cut
#Solves 4 LPs
function ab_cut_(us, vs, z, data, eps_, solver, CASE)
    const TOL = 1e-6
    m = solver==nothing ? Model() : Model(solver=solver)
    @defVar(m, a[1:size(data, 2)])
    @defVar(m, b)

    #define the RHS
    const N = size(data, 1); const d = size(data, 2)
    @defVar(m, t[1:N] >= 0)
    for i = 1:N
        @addConstraint(m, t[i] >= z*(dot(a, vec(data[i, :])) - b)/N)
    end

    #a,b in a box
    @defVar(m, abs_a[1:d])
    for i = 1:d
        @addConstraint(m, abs_a[i] >=  a[i])
        @addConstraint(m, abs_a[i] >= -a[i])
    end
    #Altered from original paper to match the new sampling technique
    @addConstraint(m, sum{abs_a[i], i=1:d} <= 1.)
 
    data_max = maximum(abs(data))

    if CASE == :1
        @addConstraint(m, dot(a, vs) - (z-1) * b >= 0)
        @addConstraint(m, dot(a, us) - b >= 0)

        #implied bounds on b to help maintain numerical stability
        @addConstraint(m, b >= -data_max)

        @setObjective(m, :Max, dot(a, vs) - (z-1)*b + dot(a, us) - b - sum{t[i], i=1:N})
    elseif CASE == :2
        @addConstraint(m, dot(a, vs) - (z-1) * b >= 0)
        @addConstraint(m, dot(a, us) - b <= 0)

        #implied bounds on b to help maintain numerical stability
        @addConstraint(m, (z-1) * b <= norm(vs, Inf))
        if abs(z-1) <= TOL
            @addConstraint(m, b <= data_max)
        end

        @setObjective(m, :Max, dot(a, vs) - (z-1)*b - sum{t[i], i=1:N})
    elseif CASE == :3
        @addConstraint(m, dot(a, vs) - (z-1) * b <= 0)
        @addConstraint(m, dot(a, us) - b >= 0)

        #implied bounds on b to help maintain numerical stability
        @addConstraint(m, (z-1) * b >= -norm(vs, Inf) )
        if abs(z-1) <= TOL
            @addConstraint(m, b >= -data_max)
        end

        @setObjective(m, :Max, dot(a, us) - b - sum{t[i], i=1:N})
    end
    status = solve(m)
    status != :Optimal && error("Case $CASE optimization failed with status $status")
    getObjectiveValue(m), getValue(a), getValue(b)
end

#returns the first violation
function ab_cut(us, vs, z, data, eps_, Gamma, TOL, solver = nothing)
    obj, astar, bstar = ab_cut_(us, vs, z, data, eps_, solver, :1)
    if obj > Gamma + TOL
        return obj, astar[:], bstar
    end

    obj, astar, bstar = ab_cut_(us, vs, z, data, eps_, solver, :2)
    if obj > Gamma + TOL
        return obj, astar[:], bstar
    end

    obj, astar, bstar = ab_cut_(us, vs, z, data, eps_, solver, :3)
    if obj > Gamma + TOL
        return obj, astar[:], bstar
    end
    return Gamma, Float64[], 0.0
end

#returns zstar,ustar
#Rebuilds the model each time bc no way to remove cnsts from JuMP
function suppFcnLCX(xs, data, eps_, Gamma, cut_sense; 
                    lbounds=Float64[], ubounds=Float64[], 
                    TOL=1e-6, MAXITER=1000, trace=false, outer_solver=nothing, abcut_solver=nothing)
    #assign bounds if the user didn't give em
    const d = size(data, 2)
    data_max = maximum(abs(data))

    m = outer_solver==nothing ? Model() : Model(solver=outer_solver)
    @defVar(m, us[1:d])
    @defVar(m, vs[1:d])
    @defVar(m, 1 <= z <= 1/eps_)
    @setObjective(m, cut_sense, dot(xs, us))

    #any user bounds
    !isempty(lbounds) && @addConstraint(m, lbounds[i] <= us[i=1:d])
    !isempty(ubounds) && @addConstraint(m, ubounds[i] >= us[i=1:d])

    #Implied bounds on u,v, added for stability
    for ix = 1:d
        @addConstraint(m,  vs[ix] <= data_max * z + Gamma)
        @addConstraint(m, -vs[ix] <= data_max * z + Gamma)
        @addConstraint(m,  us[ix] <= data_max * z + Gamma)
        @addConstraint(m, -us[ix] <= data_max * z + Gamma)
    end

    status = solve(m)
    status != :Optimal && error("LCX support init solve failed: $status")

    uvals = getValue(us)
    vvals = getValue(vs)
    zval  = getValue(z)

    obj, astar, bstar = ab_cut(uvals[:], vvals[:], zval, data, eps_, Gamma, TOL, abcut_solver)


    iter::Int64 = 0
    while obj > Gamma + TOL
        iter > MAXITER && error("LCX Supp: Max Iterations exceeded, $(obj - Gamma)")
        iter += 1

        @defVar(m, t[1:2] >= 0)

        @addConstraint(m, t[1] >= dot(astar[:], vs) - bstar * z + bstar)
        @addConstraint(m, t[2] >= dot(astar[:], us) - bstar)
        @addConstraint(m, t[1] + t[2] <= mean(max(data * astar[:] - bstar, 0))*z + Gamma)
        
        status = solve(m)
        status != :Optimal && error("LCX support Inner Solver Failed: $status")
        uvals = getValue(us)
        vvals = getValue(vs)
        zval  = getValue(z)

        obj, astar, bstar = ab_cut(uvals[:], vvals[:], zval, data, eps_, Gamma, TOL, abcut_solver)

        trace && println("Iter: $iter \t $(obj-Gamma) \t $(astar') \t $(bstar)") # \t $(uvals[:]) \t $(vvals[:]) \t $zval
    end
    getObjectiveValue(m), getValue(us[:])    
end

########################
type LCXOracle <: AbstractOracle
    data::Array{Float64, 2}
    eps_::Float64
    Gamma::Float64
    cut_tol::Float64         # 1e-6
    max_iter::Int64          # 100
    ab_cut_tol::Float64      # 1e-6
    abcut_solver::Union(AbstractMathProgSolver, Nothing)
    outer_solver::Union(AbstractMathProgSolver, Nothing)
    
    lbounds::Vector{Float64}
    ubounds::Vector{Float64}

    # Other options
    debug_printcut::Bool
    trace::Bool

end

#Preferred Interface
function LCXOracle(data, eps_, delta; 
                    cut_tol=1e-6, max_iter=100, ab_cut_tol=1e-6, 
                    lbounds=Float64[], ubounds=Float64[],
                    debug_printcut=false, trace=false, 
                    numSamples=int(1e4), numBoots=int(1e4), 
                    Gamma =-1.0, abcut_solver=nothing, outer_solver=nothing) 
    @assert (0 < eps_ < 1) "Epsilon invalid: $eps_"
    Gamma_ = Gamma < 0 ? calc_ab_thresh(data, delta, numBoots, numSamples) : Gamma
    LCXOracle(data, eps_, Gamma_, cut_tol, max_iter,
                ab_cut_tol, abcut_solver, outer_solver, lbounds, ubounds, debug_printcut, trace)
end

suppFcn(xs, w::LCXOracle, cut_sense) = 
    suppFcnLCX(xs, w.data, w.eps_, w.Gamma, cut_sense,  
        lbounds=w.lbounds, ubounds=w.ubounds, TOL=w.ab_cut_tol, MAXITER=w.max_iter, 
        trace=w.trace, outer_solver=w.outer_solver, abcut_solver=w.abcut_solver)

# JuMPeR alerting us that we're handling this contraint
registerConstraint(w::LCXOracle, rm::Model, ind::Int, prefs) = 
    ! get(prefs, :prefer_cuts, true) && error("Only cutting plane supported")

function setup(w::LCXOracle, rm::Model, prefs)
    rd = JuMPeR.getRobust(rm)
    @assert (rd.numUncs == size(w.data, 2)) "Num Uncertainties $(rd.numUncs) doesn't match columns in data $(size(w.data, 2))"
    @assert (length(rd.uncertaintyset) == 0) "Auxiliary constraints on uncertainties not yet supported"

    #set the bounds on the object now.  
    #Shallow copy since they should be "const"
    #only really want them if they're not uniformly Infinite
    if ! all(rd.uncLower .== -Inf)
        w.lbounds = rd.uncLower
    end
    if ! all(rd.uncUpper .== Inf)
        w.ubounds = rd.uncUpper
    end

    # Extract preferences we care about
    w.debug_printcut = get(prefs, :debug_printcut, false)
    w.cut_tol        = get(prefs, :cut_tol, w.cut_tol)
end

function generateCut(w::LCXOracle, m::Model, rm::Model, inds::Vector{Int}, active=false)
    new_cons = {}
    rd = JuMPeR.getRobust(rm)

    for ind in inds
        con = JuMPeR.get_uncertain_constraint(rm, ind)
        cut_sense, xs, lhs_const = JuMPeR.build_cut_objective(rm, con, m.colVal)
        d = size(w.data, 2)
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
generateReform(w::LCXOracle, m::Model, rm::Model, inds::Vector{Int}) = 0
