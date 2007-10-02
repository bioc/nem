PhiDistr <- function(Phi, Pm, a=1, b=ncol(Phi)^2){	
	d <- sum(abs(Phi - Pm)) 	# take mean instead of sum to make it numerically more stable		
	a/(2*b)*(1 + d/b)^(-a-1)
}