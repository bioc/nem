PhiDistr <- function(Phi, Pm, a=1, b=0.5){	
	d <- sum(abs(Phi - Pm)) 		
	a/(2*b)*(1 + d/b)^(-a-1)
}