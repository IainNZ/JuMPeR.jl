using JuMPeR
using Base.Test

function Test1(cuts, rule)
    println(" Test1 $cuts $rule")
    
    sale_price = 3.0

    m = RobustModel(solver=solver)
    @defVar(m, buy[1:2] >= 0)
    @defVar(m, sold[1:2] >= 0)
    @defVar(m, profit)
    @defUnc(m, demand[1:2])
    setName(buy[1],"buy[1]")
    setName(buy[2],"buy[2]")
    setName(sold[1],"sold[1]")
    setName(sold[2],"sold[2]")
    setName(demand[1], "d1")
    setName(demand[2], "d2")

    @setObjective(m, Max, profit)
    addConstraint(m, profit <= -(buy[1] + buy[2]) + sale_price*(sold[1] + sold[2]))
    addConstraint(m, sold[1] <= buy[1])
    addConstraint(m, sold[2] <= buy[2])
    addConstraint(m, sold[1] <= demand[1])
    addConstraint(m, sold[2] <= demand[2])

    if rule == :DictDict
        setAdapt!(sold, :Affine, [demand])
    elseif rule == :VarDict
        setAdapt!(sold[1], :Affine, [demand])
        setAdapt!(sold[2], :Affine, [demand])
    elseif rule == :DictVar
        setAdapt!(sold, :Affine, [demand[1], demand[2]])
    elseif rule == :VarVar
        setAdapt!(sold[1], :Affine, [demand[1], demand[2]])
        setAdapt!(sold[2], :Affine, [demand[1], demand[2]])
    end

    addConstraint(m,  (demand[1] - 20)/10 + (demand[2] - 10)/5 <= 1)
    addConstraint(m,  (demand[1] - 20)/10 - (demand[2] - 10)/5 <= 1)
    addConstraint(m, -(demand[1] - 20)/10 + (demand[2] - 10)/5 <= 1)
    addConstraint(m, -(demand[1] - 20)/10 - (demand[2] - 10)/5 <= 1)

    solveRobust(m, prefer_cuts=cuts)
    @test_approx_eq getValue(profit) 35.0
    @test_approx_eq getValue(buy[1]) 15.0
    @test_approx_eq getValue(buy[2]) 10.0
end

for cut_pref in [true, false]
    for rule in [:DictDict, :VarDict, :DictVar, :VarVar]
        Test1(cut_pref, rule)
    end
end