dd_psi = function(y){
	out_psi = 0.5 + besselI(y, 2)/(2 * besselI(y, 0)) - (besselI(y, 1)/besselI(y, 0))^2
	return(out_psi)
}

vm_lb = function(n, h, a, b){
	out = 1 - 1/sqrt(h * (n * dd_psi(b) - c0 / (a^2)))
	return(out)
}