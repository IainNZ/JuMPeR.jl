#############################################################################
# JuMPeR
# Julia for Mathematical Programming - extension for Robust Optimization
# See http://github.com/IainNZ/JuMPeR.jl
#############################################################################
# scenario.jl
# Logic for scenarios, which are samples from the uncertainty set, and can
# be both provided and retrieved at optimality.
#############################################################################

export Scenario, setUncValue, getUncValue, isBinding
type Scenario
    data::Vector{Float64}  # Using NaN as undefined
    binding::Bool
end
setUncValue(s::Scenario, u::Uncertain, v::Float64) = (s.data[u.unc] = v)
getUncValue(s::Scenario, u::Uncertain) = (s.data[u.unc])
isBinding(s::Scenario) = s.binding


# addScenario
# Provide a scenario as either a dictionary or a Scenario type
export addScenario
function addScenario(m::Model, data::Dict)
    scen = Scenario(fill(NaN,getRobust(m).numUncs), false)
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

# scen_satisfies_con
# Internal function. Given a constraint and a scenario, does the scenario
# define all the uncertainties that appear in the constraint?
function scen_satisfies_con(scen::Scenario, con::UncConstraint)
    # Variable part
    for var_ind = 1:length(con.terms.coeffs)
        coeff = con.terms.coeffs[var_ind]
        for unc_ind = 1:length(coeff.vars)
            isnan(scen.data[coeff.vars[unc_ind].unc]) && return false
        end
    end
    # Non variable part
        coeff = con.terms.constant
        for unc_ind = 1:length(coeff.vars)
            isnan(scen.data[coeff.vars[unc_ind].unc]) && return false
        end
    return true
end

# scen_to_vec [internal function]
scen_to_vec(scen::Scenario) = scen.data

# cut_to_scen [internal function]
cut_to_scen(unc_val::Vector{Float64}, binding=false) = Scenario(unc_val,binding)