moduleNetwork <- function(D,type="mLL",Pe=NULL,Pm=NULL,lambda=0,para=NULL,hyperpara=NULL,selEGenes=FALSE,verbose=TRUE){				
	Sgenes = unique(colnames(D))	
	n <- length(Sgenes)	
	P <- sapply(Sgenes, function(x) rowSums(D[, colnames(D) == x, drop=FALSE]))	
	cat("Estimating module network of",n,"S-genes (lambda =", lambda,")...\n\n")
	C = cor(P)			
	if(is.null(Pe)){
		Pe <- matrix(1/n,nrow=nrow(D),ncol=n)	
		colnames(Pe) <- Sgenes			
	}		
	modeltotal <- matrix(0,ncol=n,nrow=n)
	colnames(modeltotal) <- Sgenes
	rownames(modeltotal) <- Sgenes	
	
	Dis0 = as.dist(1-C)
# 	Dis0 = dist(t(P))						
	hc <- hclust(Dis0,method="average")			
	clusterData <- function(h){
		cl <- cutree(hc,h=h)
		T = table(cl)
		list(clust = cl, sizes = T, height=h)
	}	
	heights = sort(hc$height,decreasing=TRUE)	
	clust <- lapply(heights, clusterData)		
	maxsizes = sapply(clust, function(x) max(x$sizes))
	minsizes = sapply(clust, function(x) min(x$sizes))
	hmin = min(which(maxsizes < 5))
	hmax = min(which(minsizes < 5))	
	idx = list()
	for(h in hmax:hmin){		
		if(verbose){		
			plot(hc)	
			if(h > 1)
				rect.hclust(hc,h=clust[[h-1]]$height)							
		}									
		small = which(clust[[h]]$sizes < 5)		
		for(c in small){
			variables = which(clust[[h]]$clust == c)			
			variables = setdiff(variables,unlist(idx)) # remove those, which have already been part of a module
			if(length(variables) > 0){
				if(verbose){
					if(h > 1)
						rect.hclust(hc,h=clust[[h-1]]$height,which=c,border="green")	
				}
				idx = c(idx, list(variables))		
				vars <- colnames(P)[variables]							
				varidx <- sapply(vars,function(x) which(colnames(D) %in% x))	
				modeltotal <- moduleNetwork.aux(D[,varidx,drop=FALSE],modeltotal,variables,Pe=Pe[,variables],Pm=Pm[variables,variables],lambda=lambda,para=para,hyperpara=hyperpara,type=type,selEGenes=selEGenes,verbose=verbose)
			}	
		}		
	}
	if(length(idx) > 1){	
		res <- connectModules(D, modeltotal, idx, Pe=Pe,Pm=Pm,lambda=lambda,para=para,hyperpara=hyperpara,type=type,selEGenes=selEGenes,verbose=verbose)
		modeltotal <- transitive.closure(modeltotal, mat=TRUE,loop=TRUE)
	}else{
		diag(modeltotal) <- 0	
		ep <- score(list(modeltotal),D,type=type,para=para,Pe=Pe,Pm=NULL,lambda=0,hyperpara=hyperpara,selEGenes=selEGenes,verbose=FALSE)		
	# output
		modeltotal <- as(modeltotal,"graphNEL")	
	     	res <- list(graph=modeltotal,mLL=ep$mLL[[1]],pos=ep$pos[[1]],mappos=ep$mappos[[1]],type=ep$type,para=para,hyperpara=hyperpara,lam=lambda)	# output: data likelihood under given model!
		class(res) <- "ModuleNetwork"
		if(verbose){
			cat("log-likelihood of model = ",res$mLL,"\n")
		}
	}
	return(res)
}

moduleNetwork.aux <- function(D,modeltotal, variables,Pe=NULL,Pm=NULL,lambda=0,para=NULL,hyperpara=NULL,type="mLL",selEGenes=FALSE,verbose=TRUE){
	if(verbose){
		cat("estimating network of genes:\n")
		cat(variables,"\n")
	}
	n <- length(variables)
	if(n > 1){			
		models <- enumerate.models(n,name=unique(colnames(D)),verbose)
		sco <- score(models,D,type=type,para=para,hyperpara=hyperpara,Pe=Pe,Pm=Pm,lambda=lambda,selEGenes=selEGenes,verbose=verbose,graphClass="matrix")
		modellocal <- sco$graph
	}
	else
		modellocal <- 0
	modeltotal[variables,variables] <- modellocal
	modeltotal
}

connectModules <- function(D, Phi, modules, type="mLL",Pe=NULL,Pm=NULL,lambda=0,para=NULL,hyperpara=NULL,selEGenes=FALSE,verbose=TRUE){
	Sgenes = unique(colnames(D))
	n <- length(Sgenes)	
	if(verbose){	
		dev.off()
		cat("Connecting modules using constraint greedy hillclimbing ...\n\n")	
	}
	sco0 <- score(list(Phi),D,type=type,para=para,hyperpara=hyperpara,Pe=Pe,Pm=Pm,lambda=lambda,selEGenes=selEGenes,verbose=verbose,graphClass="matrix")$mLL			
	BetweenModules = matrix(0,n,n)
	dimnames(BetweenModules) = dimnames(Phi)
	for(i in 1:length(modules)){
		mod1 = modules[[i]]	
		if(i < length(modules)){
			for(j in (i+1):length(modules)){
				mod2 = modules[[j]]				
				BetweenModules[mod1,mod2] = 1
				BetweenModules[mod2,mod1] = 1
			}
		}
	}	
	finished <- FALSE
	while(!finished){		
# 		propose new edges between modules									
		idx = which(BetweenModules == 1)
		if(length(idx) > 0){	
			models <- list()
			for(i in 1:length(idx)){ # test all possible new edges
				Phinew = Phi
				Phinew[idx[i]] = 1
				Phinew = transitive.closure(Phinew, mat=TRUE, loops=TRUE)	
				models[[i]] <- Phinew				
			}				
			models <- unique(models)
			sconew <- score(models,D,type=type,para=para,hyperpara=hyperpara,Pe=Pe,Pm=Pm,lambda=lambda,selEGenes=selEGenes,verbose=verbose,graphClass="matrix")			
			maidx = which.max(sconew$mLL)
			if(sconew$mLL[maidx] > sco0){
				newedges = which(sconew$graph - Phi == 1)
				BetweenModules[newedges] = 0
				sco0 <- sconew$mLL[maidx]
				Phi <- sconew$graph
			}
			else # otherwise no improving edge could be inserted
				finished <- TRUE
		}
		else
			finished <- TRUE	
	}
	ep <- score(list(Phi),D,type=type,para=para,Pe=Pe,Pm=NULL,lambda=0,hyperpara=hyperpara,selEGenes=selEGenes,verbose=FALSE,graphClass="matrix")
	diag(Phi) <- 0
	Phi <- as(Phi,"graphNEL")	
    	res <- list(graph=Phi,mLL=ep$mLL[[1]],pos=ep$pos[[1]],mappos=ep$mappos[[1]],type=ep$type,para=para,hyperpara=hyperpara,lam=lambda)	# output: data likelihood under given model!
	class(res) <- "ModuleNetwork"
	if(verbose)
		cat("log-likelihood of model = ",res$mLL,"\n")
	return(res)
}