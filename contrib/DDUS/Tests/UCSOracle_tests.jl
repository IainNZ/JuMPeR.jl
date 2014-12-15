###
# UCS Oracle Tests
###
facts("suppFcnTest UCS") do
	muhat  = [0.0, 1.2]
	Gamma1 = .002
	Gamma2 = .05
	covbar = [ 1.  0.2 
	           0.2 1.3]
	
	covhat = covbar - Gamma2 * eye(2)  #covbar = covhat + Gamma2 I
	w = UCSOracle(muhat, covhat, Gamma1, Gamma2, .1)

	zstar, ustar = suppFcn([1, 1], w, :Max)
	@fact zstar => roughly(6.132331444671241)
	@fact ustar[1] => roughly(2.4661657223356204)
	@fact ustar[2] => roughly(3.66616572233562)

	zstar, ustar = suppFcn([1, 1], w, :Min)
	@fact zstar => roughly(-3.7323314446712406)
	@fact ustar[1] => roughly(-2.4661657223356204)
	@fact ustar[2] => roughly(-1.266165722335620)
end

#doesn't use any bounds
facts("portTest UCS No bounds") do
	srand(8675309); data = randn(500, 2)
	w = UCSOracle(data, .1, .1, .1)
	portTest(w, -2.4630938710200345, [0.5253456488147422, 0.4746543511852578])
end

facts("portTest2  UCS bounds") do
	srand(8675309); data = randn(500, 2)
	w = UCSOracle(data, .1, .1, .1)
	portTest(w, -2.4604223389522866, [0.530985977658841, 0.46901402234115896], TOL=1e-6,
			unc_lower=[-1e6, -1e6], unc_upper=[1e6, 1e6])
end
