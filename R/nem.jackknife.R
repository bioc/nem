nem.jackknife <- function(D, inference="nem.greedy",models=NULL,type="mLL",para=NULL,hyperpara=NULL,Pe=NULL,Pm=NULL,Pmlocal=NULL,local.prior.size=length(unique(colnames(D))),local.prior.bias=1,triples.thrsh=0.5,lambda=0,delta=1,selEGenes=FALSE,verbose=TRUE){
	inferNetwork <- function(sgenesk){	
		if(verbose)
			cat("S-genes used: ", sgenesk,"\n")
		Dtmp = D[,sapply(sgenesk,function(x) which(colnames(D) %in% x))]
		Petmp = Pe[,sgenesk]	
		Pmtmp = Pm[sgenesk,sgenesk]
		if(!is.null(Pm) & length(lambda) > 1)
			res = as(nemModelSelection(lambda,Dtmp,inference,models,type,para,hyperpara,Petmp,Pmtmp,Pmlocal,local.prior.size-1,local.prior.bias,triples.thrsh,delta,selEGenes,verbose)$graph,"matrix")
		else
			res = as(nem(Dtmp,inference,models,type,para,hyperpara,Petmp,Pmtmp,Pmlocal,local.prior.size-1,local.prior.bias,triples.thrsh,lambda,delta,selEGenes,verbose)$graph,"matrix")
		overlapJack[sgenesk,sgenesk] = overlapJack[sgenesk,sgenesk] + res
		overlapJack
	}
	Sgenes = unique(colnames(D))
	n = length(Sgenes)		
	overlapJack = matrix(0,ncol=n,nrow=n)	
	colnames(overlapJack) = Sgenes
	rownames(overlapJack) = Sgenes	
	for(i in 1:n)
		overlapJack = inferNetwork(Sgenes[-i])
	overlapJack = overlapJack / (n-1)	
	overlapJack = round(overlapJack,digits=2)
	overlapJack	
}
