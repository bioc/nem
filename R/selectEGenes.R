filterEGenes = function(Porig, D, ntop=100){
	ntop = min(nrow(Porig),ntop)
	n = ncol(D)
	I1 = apply(Porig,2,function(x) order(x)[1:ntop])
	I1 = unique(as.vector(I1))  # Länge ist abhängig von P!	
	print(paste("selecting top",ntop," genes from each list -->",length(I1),"genes total"))	
	print("Clustering E genes:")
	DMat = as.dist(1-cor(t(D[I1,])))
	hc = hclust(DMat,method="ward")
	clusterData <- function(k){		
		cl = cutree(hc,k)
		s <- summary(silhouette(cl,DMat))$avg.width
		if(is.null(s))
			s = 0
		cat(".")
		return(list(sil=s, clust=cl, k=k))
	}	
	clust = lapply(unique(pmin(n*(8:20),length(I1)-1)),clusterData)	
	maidx <- which.max(sapply(clust, function(x) x$sil))	
	kbest <- clust[[maidx]]$k
	sil <- clust[[maidx]]$sil
	clust <- clust[[maidx]]$clust						
	cat("Using ", kbest, "clusters (silhouette index:", sil ,").\n")
	I = sapply(1:kbest, function(k){
		cl = which(clust == k)	
		if(length(cl) > 2)		
			cl[order(rowSums(D[cl,]),decreasing=TRUE)[1:min(1,sum(clust==k))]]
	})
	I = unique(unlist(I))
	I = I[which(apply(D[I,],1,var) > 0)]
	cat("---> ", length(I), "E-genes selected\n")
	I
}

selectEGenes <- function(Phi,D1,D0=NULL,para=NULL,hyperpara=NULL,Pe=NULL,Pm=NULL,lambda=0,type="mLL", nEgenes=min(5*ncol(Phi), nrow(D1))){		
	mLLWrapper <- function(egene, D1, D0=NULL,a=0.05,b=0.15,Pe=NULL,type="mLL"){	
		idx <- which(rownames(D1) == egene)[1]									
		if(!is.null(Pe))	
			Pesel <- t(as.matrix(Pe[idx,]))
		else
			Pesel <- NULL		
		if(type == "FULLmLL")
			return(FULLmLL(Phi,t(as.matrix(D1[idx,])),t(as.matrix(D0[idx,])),a0=hyperpara[1],b0=hyperpara[2],a1=hyperpara[3],b1=hyperpara[4],Pe=Pesel,Pm=NULL,lambda=0)$mLL)
		else	
			return(mLL(Phi,t(as.matrix(D1[idx,])),t(as.matrix(D0[idx,])),a=para[1],b=para[2],Pe=Pesel,Pm=NULL,lambda=0,type=type)$mLL)
	}				
	if(is.null(rownames(D1)))
		rownames(D1) <- as.character(1:nrow(D1))
	L <- sapply(rownames(D1),mLLWrapper,D1=D1,D0=D0,a=a,b=b,Pe=Pe,type=type)		
	if(type %in% c("CONTmLLDens"))
		sel <- order(L,decreasing=TRUE)[1:length(which(L>0))]		
	else	
		sel <- order(L,decreasing=TRUE)[1:nEgenes]			
	if(type == "FULLmLL")
		return(c(FULLmLL(Phi,D1[sel,],D0[sel,],a0=hyperpara[1],b0=hyperpara[2],a1=hyperpara[3],b1=hyperpara[4],Pe=Pe[sel,],Pm=Pm,lambda=lambda),loglik=L))
	else
		return(c(mLL(Phi,D1[sel,],D0[sel,],a=para[1],b=para[2],Pe=Pe[sel,],Pm=Pm,lambda=lambda,type=type),loglik=L))
}

getRelevantEGenes <- function(Phi, D, nEgenes=min(5*ncol(Phi), nrow(D1)), type="mLL", para=NULL, hyperpara=NULL, Pe=NULL, Pm=NULL, lambda=0){
	# Which Sgenes were silenced?
	Sgenes <- unique(colnames(D))
	nrS <- length(Sgenes)
	
	# check that all models have S-genes as names
	fkt <- function(x,s){
		ss <- sort(s)
		c1 <- all(sort(colnames(x))==ss)
		c2 <- all(sort(rownames(x))==ss)
		return(c1 & c2)
	}	
	nrS <- length(Sgenes)    
	# if no prior is supplied:
	# assume uniform prior over E-gene positions
	if (is.null(Pe)){ 
		Pe <- matrix(1/nrS,nrow=nrow(D),ncol=nrS)
		colnames(Pe) <- Sgenes  
	}      
	if(is.null(Pm)) lambda <- 0  
	# make probability/density matrices D0 and D1  
	# nrow=#E-genes and ncol=#S-genes        
	if(type %in% c("CONTmLL","CONTmLLDens")){   	  	
		# D1[i,j] = probability/density of EFFECT at E_i when S_j was silenced  	  	  	
		D1 <- sapply(Sgenes, function(x) apply(D[, colnames(D) == x, drop=FALSE],1,mean))		
		D0 <- NULL		
	}  
	else{  	
		# D0[i,j] = how often there is NO EFFECT at E_i when S_j was silenced
		# D1[i,j] = how often there is    EFFECT at E_i when S_j was silenced  
		D0  <- matrix(0,ncol=nrS,nrow=nrow(D),dimnames=list(rownames(D),Sgenes))	
		D1  <- D0
		for (i in 1:nrS) {
			Di     <- D[,colnames(D) == Sgenes[i],drop=FALSE]		
			D0[,i] <- rowSums(Di==0)
			D1[,i] <- rowSums(Di==1) 
		}		
	}
	selectEGenes(Phi, D1, D0, para, hyperpara, Pe, Pm, lambda, type, nEgenes)
}