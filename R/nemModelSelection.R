nemModelSelection <- function(lambdas,D,inference="nem.greedy",models=NULL,type="mLL",para=NULL,hyperpara=NULL,Pe=NULL,Pmlocal=NULL,Pm=NULL,local.prior.size=length(unique(colnames(D))),local.prior.bias=1,triples.thrsh=0.5,delta=1,selEGenes=FALSE, verbose=TRUE,...){

	infer <- function(lam){					
		res <- nem(D,inference=inference,models=models,type=type,para=para,hyperpara=hyperpara,Pe=Pe,Pmlocal=Pmlocal,Pm=Pm,local.prior.size=local.prior.size,local.prior.bias=local.prior.bias,lambda=lam,delta=delta,selEGenes=selEGenes, verbose=verbose) # ACHTUNG: nem spuckt immer den MAP score aus!!!					
		if(length(res$selected) <= 1)
			res$mLL = -Inf
		else{	
			resmLL <- score(list(as(res$graph,"matrix")),D[res$selected,,drop=FALSE],type=type,para=para,Pe=Pe,Pm=NULL,lambda=0,delta=delta,hyperpara=hyperpara,verbose=FALSE)$mLL # get true mLL
			res$mLL = resmLL		
		}
		return(res)
	}	
	results <- lapply(lambdas,infer)			
	AICs <- lapply(results,network.AIC,verbose=verbose,Pm,k=log(ncol(Pm)),...)		
	winner <- results[[which.min(AICs)]]	
	if(verbose)
		cat(paste("====> chosen best model with lambda =",winner$lam,"\n"))		
	return(winner)
}
