#### 
# LCX Oracle tests
####

facts("suppFcnTest LCX") do
	srand(8675309); data = randn(100, 2)

	#VG There is a disturbing dependency on the solver here!!!!
	w = LCXOracle(data, .1, 0; Gamma =.35) 
	zstar, ustar = suppFcn([1, 1], w, :Min)
	@fact zstar => roughly(-3.1455496412236243)
	@fact ustar[1] => roughly(-1.313810356077738)
	@fact ustar[2] => roughly(-1.8317392851458862)

	zstar, ustar = suppFcn([1, 1], w, :Max)
	@fact zstar => roughly(3.447833040263504)
	@fact ustar[1] => roughly(1.5387299947906286)
	@fact ustar[2] => roughly(1.9091030454728752)

	zstar, ustar = suppFcn([-1, 1], w, :Min)
	@fact zstar => roughly(-2.8182398717617616)
	@fact ustar[1] => roughly(1.3068007701076567)
	@fact ustar[2] => roughly(-1.5114391016541049)

	zstar, ustar = suppFcn([-1, 1], w, :Max)
	@fact zstar => roughly(3.157077158917673)
	@fact ustar[1] => roughly(-1.7212256391067058)
	@fact ustar[2] => roughly(1.4358515198109671)
end


facts("portTest LCX") do
	srand(8675309); data = randn(500, 2)
	w = LCXOracle(data, .1, .2, Gamma=.2, ab_cut_tol=5e-6)
	portTest(w, -1.4293015271256384, [0.5424436347803758, 0.4575563652196242])
end