nemModelSelection <- function(lambdas,D,inference="pairwise",models=NULL,type="mLL",para=NULL,hyperpara=NULL,Pe=NULL,Pmlocal=NULL,Pm=NULL,local.prior.size=length(unique(colnames(D))),local.prior.bias=1,triples.thrsh=0.5,selEGenes=FALSE,verbose=TRUE,...){

	infer <- function(lam,D,inference,models,type,para,hyperpara,Pe,Pm,Pmlocal,local.prior.size,local.prior.bias,triples.thrsh,selEGenes,verbose){					
		res <- nem(D,inference=inference,models=models,type=type,para=para,hyperpara=hyperpara,Pe=Pe,Pmlocal=Pmlocal,Pm=Pm,local.prior.size=local.prior.size,local.prior.bias=local.prior.bias,lambda=lam,selEGenes=selEGenes,verbose=verbose)
		if(inference == "search"){
			res <- score(list(as(res$graph,"matrix")),D,type=type,para=para,Pe=Pe,Pm=NULL,lambda=0,hyperpara=hyperpara,selEGenes=selEGenes,verbose=FALSE)		
			return(res)
		}			
		return(res)
	}	
	results <- lapply(lambdas,infer,D,inference,models,type,para,hyperpara,Pe,Pm,Pmlocal,local.prior.size,local.prior.bias,triples.thrsh,selEGenes,verbose)
	AICs <- lapply(results,network.AIC,verbose=verbose,...)	
	winner <- results[[which.min(AICs)]]	
	if(verbose)
		cat(paste("====> chosen best model with lambda =",winner$lam,"\n"))
	return(winner)
}