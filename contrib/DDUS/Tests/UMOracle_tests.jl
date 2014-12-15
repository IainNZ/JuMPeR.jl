#### 
# UM Oracle tests
####

facts("suppFcnTest UM") do
	lquants = [-.05, -.25]
	uquants = [.1, .2]

	w = UMOracle(lquants, uquants, 1e-6, false)
	zstar, ustar = suppFcn([1, 1], w, :Min)
	@fact zstar =>  roughly(-0.3)
	@fact ustar[1] => roughly(-0.05)
	@fact ustar[2] => roughly(-0.25)

	zstar, ustar = suppFcn([1, 1], w, :Max)
	@fact zstar => roughly(.3)
	@fact ustar[1] => roughly(.1)
	@fact ustar[2] => roughly(.2)

	zstar, ustar = suppFcn([-1,  1], w, :Min)
	@fact zstar => roughly(-.35)
	@fact ustar[1] => roughly(.1)
	@fact ustar[2] => roughly(-.25)

	zstar, ustar = suppFcn([-1,  1], w, :Max)
	@fact zstar => roughly(.25)
	@fact ustar[1] => roughly(-.05)
	@fact ustar[2] => roughly(.2)
end

facts("portTest for UM") do
	srand(8675309); data = randn(500, 2)
	w = UMOracle(data, -3.1*ones(2), 2.9 * ones(2), .1, .2)
	portTest(w, -1.7725987075971736, [1, 0])
end
