nem.greedyMAP <- function(D, Pe=NULL, Pm=NULL, lambda=0, delta=1, tol=1e-4, trans.close=TRUE, verbose=TRUE){		
	Sgenes = unique(colnames(D))
	n <- length(Sgenes)		
	cat("Alternating optimization for",n,"S-genes (lambda =", lambda,")...\n\n")
	Theta = apply(D,1, function(e) (e == max(e)) & (e>0))*1		
	best = NULL
	converged = FALSE
	i = 1	
	while(!converged){		
		Phi = apply(Theta%*%D, 1, function(e) (e > 0))*1							
		sco <- score(list(Phi),D,type="CONTmLLMAP",Pe=Pe,Pm=Pm,lambda=lambda,delta=delta,verbose=verbose)
		if(verbose){
			cat("iteration ", i, ": likelihood = ", sco$mLL, "\n")
		}
		if(!is.null(best))		
			converged = (sco$mLL/sco0$mLL - 1 < tol)
		if(converged)
			break;		
		Theta = apply(sco$pos[[1]],1,function(e) (e == max(e)) & (e>0))*1	
		Theta = Theta[1:(nrow(Theta)-1),]
		sco0 = sco			
		if(is.null(best) || (sco$mLL > best$mLL))
			best = sco
		i = i + 1
	}			
	sco = best
	Phi = as(sco$graph, "matrix")	
	if(trans.close)	
		Phi = closest.transitive.greedy(Phi, verbose)
	ep <- score(list(Phi),D,type="CONTmLLMAP",Pe=Pe,Pm=Pm,lambda=lambda,delta=delta,verbose=verbose)
    	res <- list(graph=ep$graph,mLL=ep$mLL[[1]],pos=ep$pos[[1]],mappos=ep$mappos[[1]],type=ep$type[[1]],para=ep$para,hyperpara=ep$hyperpara,lam=lambda,selected=ep$selected,delta=delta,LLperGene=ep$LLperGene[[1]])	# output: data likelihood under given model!	
	class(res) <- "nem.greedyMAP"	
	return(res)
}
