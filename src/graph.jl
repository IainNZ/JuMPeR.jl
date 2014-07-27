#############################################################################
# JuMPeR
# Julia for Mathematical Programming - extension for Robust Optimization
# See http://github.com/IainNZ/JuMPeR.jl
#############################################################################
# graph.jl
# Graph-related algorithms to more efficiently organize the optimization
#############################################################################

function detect_components(num_unc, uncset)
    unc_to_comp = zeros(Int, num_unc)
    con_to_comp = zeros(Int, length(uncset))

    function replace_component(from, to)
        # Change all `from` to `to`
        for k in 1:num_unc
            unc_to_comp[k] != from && continue
            unc_to_comp[k] = to
        end
        for i in 1:length(uncset)
            con_to_comp[i] != from && continue
            con_to_comp[i] = to
        end
    end


    num_components = 0
    for con_ind in 1:length(uncset)
        lhs = uncset[con_ind].terms
        length(lhs.vars) == 0 && continue
        new_component_needed = true
        row_component = 0
        for j = 1:length(lhs.vars)
            unc      = lhs.vars[j].unc
            unc_comp = unc_to_comp[unc]
            unc_comp == 0 && continue
            # Something in this row has a component already - convert
            # everything in this row to that component.
            for k = 1:length(lhs.vars)
                k_unc      = lhs.vars[k].unc
                k_unc_comp = unc_to_comp[k_unc]
                if k_unc_comp == 0
                    unc_to_comp[k_unc] = unc_comp
                elseif k_unc_comp != unc_comp
                    replace_component(k_unc_comp, unc_comp)
                end
            end
            con_to_comp[con_ind] = unc_comp
            new_component_needed = false
            break
        end
        if new_component_needed
            # Entire row was component-less, create new component
            num_components += 1
            for v in lhs.vars
                unc_to_comp[v.unc] = num_components
            end
            con_to_comp[con_ind] = num_components
        end
    end

    # Compress components down
    d = [v => u for (u,v) in enumerate(sort(collect(Set(unc_to_comp))))]

    return Int[d[i] for i in unc_to_comp], Int[d[i] for i in con_to_comp]
end