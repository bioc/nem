nem.bootstrap <- function(D,nboot=1000,inference="nem.greedy",models=NULL,type="mLL",para=NULL,hyperpara=NULL,Pe=NULL,Pm=NULL,Pmlocal=NULL,local.prior.size=length(unique(colnames(D))),local.prior.bias=1,triples.thrsh=0.5,lambda=0,delta=1,selEGenes=FALSE,verbose=TRUE){
		
	inferNetwork <- function(boot){
		Dtmp = D[boot,]						
		PeBoot = Pe[boot,]		
		if(!is.null(Pm) & length(lambda) > 1)
			res = as.vector(as(nemModelSelection(lambda,Dtmp,inference,models,type,para,hyperpara,PeBoot,Pm,Pmlocal,local.prior.size,local.prior.bias,triples.thrsh,delta,selEGenes,verbose)$graph,"matrix"))
		else
			res = as.vector(as(nem(Dtmp,inference,models,type,para,hyperpara,PeBoot,Pm,Pmlocal,local.prior.size,local.prior.bias,triples.thrsh,lambda,delta,selEGenes,verbose)$graph,"matrix"))
		res
	}
	results = bootstrap(1:nrow(D),nboot,theta=inferNetwork)$thetastar
	Sgenes = unique(colnames(D))
	n = length(Sgenes)
	overlapBoot = rowMeans(results)
	overlapBoot = matrix(round(overlapBoot,digits=2),ncol=n,nrow=n)
	colnames(overlapBoot) = Sgenes
	rownames(overlapBoot) = Sgenes
	overlapBoot
}
