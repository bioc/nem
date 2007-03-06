# 	get AIC for network
network.AIC = function(network, k=2, verbose=TRUE){	
	M <- as(network$graph,"matrix")
	M <- M - diag(diag(M))
	npar <- sum(M != 0)		
	AIC <- -2*network$mLL + k*npar
	if(verbose) cat(paste("==> AIC ( lambda = ",network$lam,") = ",AIC,"( #param =", npar,")===============\n"))
	return(AIC)
}
