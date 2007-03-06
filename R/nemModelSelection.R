nemModelSelection <- function(lambdas,D,inference="pairwise",models=NULL,type="mLL",para=NULL,hyperpara=NULL,Pe=NULL,Pmlocal=NULL,Pm=NULL,local.prior.size=length(unique(colnames(D))),local.prior.bias=1,triples.thrsh=0.5,verbose=TRUE,...){

	infer <- function(lam,D,inference,models,type,para,hyperpara,Pe,Pm,Pmlocal,local.prior.size,local.prior.bias,triples.thrsh,verbose){			
		res <- nem(D,inference=inference,models=models,type=type,para=para,hyperpara=hyperpara,Pe=Pe,Pmlocal=Pmlocal,Pm=Pm,local.prior.size=local.prior.size,local.prior.bias=local.prior.bias,lambda=lam,verbose=verbose)
		if(inference == "search"){
			s <- unlist(res$mLL)
			ep <- res$pos
			map <- res$mappos
			best <- which.max(s)			
			res<- list(graph=res$graph, mLL=s[[best]], pos=ep[[best]], mappos=map[[best]], type=type, para=para, hyperpara=hyperpara, lam=lambda)
			class(res) <- "score"
		}	
		return(res)
	}	
	results <- lapply(lambdas,infer,D,inference,models,type,para,hyperpara,Pe,Pm,Pmlocal,local.prior.size,local.prior.bias,triples.thrsh,verbose)
	AICs <- lapply(results,network.AIC,...)	
	winner <- results[[which.min(AICs)]]	
	if(verbose)
		cat(paste("====> chosen best model with lambda =",winner$lam,"\n"))
	return(winner)
}