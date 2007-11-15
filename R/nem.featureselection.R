nem.featureselection <- function(D,inference="nem.greedy",models=NULL,type="mLL",para=NULL,hyperpara=NULL,Pe=NULL,Pm=NULL,Pmlocal=NULL,local.prior.size=length(unique(colnames(D))),local.prior.bias=1,triples.thrsh=0.5,lambda=0,verbose=TRUE, tol=1e-4){
	
	if(type == "CONTmLLRatio"){
		Sgenes = unique(colnames(D))
		nrS = length(Sgenes)
		 if (is.null(Pe)){ 	
			Pe <- matrix(1/nrS,nrow=nrow(D),ncol=nrS)
			colnames(Pe) <- Sgenes  		
  		}    		
		deltaseq = 1:10
		results = sapply(deltaseq, function(d){			
			net = nem(D, inference, models, type, para, hyperpara, Pe, Pm, Pmlocal, local.prior.size, local.prior.bias, triples.thrsh, lambda, delta=d, selEGenes=FALSE,verbose=verbose)	
			if(verbose)
				cat("Selected E-genes (delta = ", d, ", AIC = ", -net$mLL + length(net$selected), "):", sort(net$selected),"\n")				
			net
		})		
 		s = -unlist(results["mLL",]) + sapply(results["selected",], length)
 		winner = which.min(s)		
		net = results[,winner]		
	}
	else{
		converged = FALSE
		while(!converged){
			net = nem(D, inference, models, type, para, hyperpara, Pe, Pm, Pmlocal, local.prior.size, local.prior.bias, triples.thrsh, lambda, selEGenes=FALSE,verbose=verbose)
			sel = getRelevantEGenes(as(net$graph,"matrix"), D, para, hyperpara, Pe, Pm, lambda, delta=1e-10, type=type)
			if(verbose)
				cat("Old Likelihood = ", net$mLL, "new likelihood = ", sel$mLL, "\nSelected E-genes:", sort(sel$selected),"\n")			
			converged = (abs(net$mLL - sel$mLL)/abs(net$mLL) < tol)
			if(converged)
				break
			D = D[sel$selected,]
		}
		net$selected = sel$selected
		net$mLL = sel$mLL		
	}
	if(verbose)
		cat("Converged! ", length(net$selected),"selected E-genes\n")
	net
}
