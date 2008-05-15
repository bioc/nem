nem.jackknife <- function(D, thresh=0.5, inference="nem.greedy",models=NULL,type="mLL",para=NULL,hyperpara=NULL,Pe=NULL,Pm=NULL,Pmlocal=NULL,local.prior.size=length(unique(colnames(D))),local.prior.bias=1,triples.thrsh=0.5,lambda=0,delta=1,selEGenes=FALSE,verbose=TRUE){
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
	overlapJack = overlapJack / (n-2)	
	overlapJack = round(overlapJack,digits=2)			
	res = nem(D,models=list((overlapJack>thresh)*1),inference="search",type=type,para=para,Pe=Pe,Pm=Pm,lambda=lambda,delta=delta,hyperpara=hyperpara,selEGenes=selEGenes, verbose=verbose)	
	res$pos = res$pos[[1]]
	res$mappos = res$mappos[[1]]
	res$mLL = res$mLL[[1]]
	res$LLperGene = res$LLperGene[[1]]
	g = res$graph
	edgeDataDefaults(g, "label") = 1	
	edgeDataDefaults(g, "weight") = 1
	for(s1 in Sgenes){
		for(s2 in Sgenes){
			if(s2 %in% unlist(adj(g, s1))){
				edgeData(g, from = s1, to = s2, attr = "weight") = overlapJack[s1,s2]			
				edgeData(g, from = s1, to = s2, attr = "label") = overlapJack[s1,s2]
			}
		}
	}
	res$graph = g
	class(res) <- "nem.jackknife"
	res
}
