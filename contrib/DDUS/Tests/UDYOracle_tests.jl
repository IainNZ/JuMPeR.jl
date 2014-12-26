###
# UDY Oracle Tests
###
#doesn't use any bounds
facts("portTest UDY") do
	srand(8675309); data = randn(500, 2)
	w = UDYOracle(data, .1, .2)
	portTest(w, -2.450295808139833, [0.5372084714704519, 0.4627915285295481], TOL=1e-6,
			unc_lower=[-10, -10], unc_upper=[10, 10])
end

