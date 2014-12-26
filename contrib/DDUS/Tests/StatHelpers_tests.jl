###
# Helper Tests
###
include("../helpers.jl")

facts("bootstrapping tests") do
	#single variate data
	srand(8765309)
	data = randn(1000)
	@fact boot(data, mean, .9, 10000) =>roughly(0.10809317304982774, 1e-8)

	#multivariate data
	srand(8765309)
	data = randn(1000, 3)
	f= x->mean(minimum(x, 1)) # a nonsense function
	@fact boot(data, f, .9, 10000) => roughly(-2.8777884239388682, 1e-8)
end

facts("t-approximations") do
	srand(8765309)
	data = randn(1000, 3)
	out = calcMeansT(data, .1)
	@fact out[1] => roughly(0.004451891728022329, 1e-8)
	@fact out[2] => roughly(0.06592683689569365, 1e-8)

	out = calcMeansT(data, .1, joint=false)
	@fact out[1] => roughly(0.011243060792593858, 1e-8)
	@fact out[2] => roughly(0.05913566783112214, 1e-8)
end

facts("LogMeanTest") do
	srand(8675309); data = rand(10)
	@fact logMeanExp(2., data) => roughly(1.1566289640224618, 1e-8)
	@fact logMeanExp(-1., data) => roughly(-0.4493802878492036, 1e-8)
end

facts("fwdBackSigsTest") do
	srand(8675309); data = randn(100)
	sigfwd, sigback = calcSigsBoot(data, .1, 10000)
	@fact sigfwd => roughly(1.0918195335954186, 1e-10)
	@fact sigback => roughly(1.08588753435207, 1e-10)
end

facts("KSGammaTest") do
	@fact KSGamma(.1, 1000) => roughly(0.0385517413380297)
end

facts("boot_mu_test") do
	srand(8675309); data = randn(100)
	@fact boot_mu(data, .1, 100) => roughly(0.15990283885940632)
end

facts("boot_sigma_test") do
	srand(8675309); data = randn(100)
	@fact boot_sigma(data, .1, 100) => roughly(0.196348437082907)
end

facts("ab_thresh_test") do
	srand(8675309); data = randn(100, 3)
	@fact calc_ab_thresh(data, .2, 100, 100) => roughly(0.1777841615666656)
end


