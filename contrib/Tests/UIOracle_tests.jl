#### 
# UI Oracle tests
####
facts("suppFcnTest UI") do
	srand(8675309); data = rand(200, 2)
	w = UIOracle(data, [0, 0], [1, 1], .1, .2) 

	# zstar, ustar = suppFcnUI([1, 1], data, [0, 0], [1, 1], log(1/.1), .001)
	zstar, ustar = suppFcn([1, 1], w, :Max)	
	@fact zstar => roughly(1.7738298026385046, 1e-8)
	@fact ustar[1] => roughly(0.8885733427479442, 1e-8)
	@fact ustar[2] => roughly(0.8852564598905605, 1e-8)

	# zstar, ustar = suppFcnUI([1, -1], data, [0, 0], [1, 1], log(1/.1), .001)
	zstar, ustar = suppFcn([1, -1], w, :Max)
	@fact zstar => roughly(0.7921571113907224, 1e-8)
	@fact ustar[1] => roughly(0.8976922795298192, 1e-8)
	@fact ustar[2] => roughly(0.10553516813909677, 1e-8)

	# zstar, ustar = suppFcnUI([-1, 1], data, [0, 0], [1, 1], log(1/.1), .001)
	zstar, ustar = suppFcn([-1, 1], w, :Max)
	@fact zstar => roughly(0.7753270721739416, 1e-8)
	@fact ustar[1] => roughly(0.11139129191787517, 1e-8)
	@fact ustar[2] => roughly(0.8867183640918167, 1e-8)
end

facts("portTest UI") do
	srand(8675309)
	data = rand(500, 2)
	w = UIOracle(data, [0., 0.], [1., 1.], .1, .2) 
	portTest(w, 0.12335160142638069, [0.5639829612402062, 0.4360170387597938])
end
