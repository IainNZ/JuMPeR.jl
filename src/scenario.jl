#############################################################################
# JuMPeR
# Julia for Mathematical Programming - extension for Robust Optimization
# See http://github.com/IainNZ/JuMPeR.jl
#############################################################################
# scenario.jl
# Logic for scenarios, which are samples from the uncertainty set, and can
# be both provided and retrieved at optimality.
#############################################################################

export Scenario, setUncValue, getUncValue
type Scenario
    data::Vector{Float64} #::Dict{Uncertain,Float64}
end
setUncValue(s::Scenario, u::Uncertain, v::Float64) = (s.data[u.unc] = v)
getUncValue(s::Scenario, u::Uncertain) = (s.data[u.unc])


# addScenario
# Provide a scenario as either a dictionary or a Scenario type
export addScenario
function addScenario(m::Model, data::Dict{Uncertain,Float64})
    scen = Scenario(zeros(getRobust(m).numUncs))
    for u in keys(data)
        scen.data[u.unc] = data[u]
    end
    addScenario(m, scen)
end
addScenario(m::Model, scen::Scenario) = push!(getRobust(m).scenarios, scen)


# getScenario
# Given a constraint reference, get a Scenario
export getScenario
function getScenario(cr::ConstraintRef{UncConstraint})
    ac = getRobust(cr.m).activecuts[cr.idx]
    length(ac) == 0 && return nothing
    return ac[1]
end

#############################################################################

# scen_to_vec
# Internal function. Take scenario, turn it into a numUncs vector
# with values filled out
function scen_to_vec(rm::Model, scen::Scenario)
    #=unc_val = zeros(getRobust(rm).numUncs)
    for unc_obj in keys(scen.data)
        unc_val[unc_obj.unc] = getUncValue(scen, unc_obj)
    end
    return unc_val=#
    return scen.data
end

# vec_to_scen
# Internal function. Take a constraint and an associated vector
# of uncertainty values, and turn it a Scenario object
function vec_to_scen(unc_con::UncConstraint,
                     unc_val::Vector{Float64})
    #data = Dict{Int,Float64}()
    #= rm = nothing

    # Variable part
    for var_ind = 1:length(unc_con.terms.vars)
        coeff = unc_con.terms.coeffs[var_ind]
        for unc_ind = 1:length(coeff.vars)
            rm = coeff.vars[unc_ind].m
            #unc = coeff.vars[unc_ind].unc
            #data[unc] = get(data, unc, 0.0) + unc_val[unc]
        end
    end
    # Non variable part
        coeff = unc_con.terms.constant
        for unc_ind = 1:length(coeff.vars)
            rm = coeff.vars[unc_ind].m
            #unc = coeff.vars[unc_ind].unc
            #data[unc] = get(data, unc, 0.0) + unc_val[unc]
        end

    #return Scenario([Uncertain(rm,unc) => data[unc] for unc in keys(data)])

    return Scenario([Uncertain(rm,unc) => unc_val[unc] for unc=1:length(unc_val)])=#
    return Scenario(unc_val)
end