nem.consensus <- function(D,thresh=0.5, nboot=1000,inference="nem.greedy",models=NULL,type="mLL",para=NULL,hyperpara=NULL,Pe=NULL,Pm=NULL,Pmlocal=NULL,local.prior.size=length(unique(colnames(D))),local.prior.bias=1,triples.thrsh=0.5,lambda=0,delta=1,selEGenes=FALSE,verbose=TRUE){
	overlapBoot = as(nem.bootstrap(D,thresh, nboot,inference,models,type,para,hyperpara,Pe,Pm,Pmlocal,local.prior.size,local.prior.bias,triples.thrsh,lambda,delta,selEGenes,verbose)$graph,"matrix")
	overlapJack = as(nem.jackknife(D,thresh, inference,models,type,para,hyperpara,Pe,Pm,Pmlocal,local.prior.size,local.prior.bias,triples.thrsh,lambda,delta,selEGenes,verbose)$graph,"matrix")
	consens = ((overlapBoot > thresh) & (overlapJack > thresh))*1
	res = nem(D,models=list(consens),inference="search",type=type,para=para,Pe=Pe,Pm=Pm,lambda=lambda,delta=delta,hyperpara=hyperpara,selEGenes=selEGenes, verbose=verbose)
	res$pos = res$pos[[1]]
	res$mappos = res$mappos[[1]]
	res$mLL = res$mLL[[1]]
	res$LLperGene = res$LLperGene[[1]]
	class(res) <- "nem.consensus"
	res
}
