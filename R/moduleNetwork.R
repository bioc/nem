moduleNetwork <- function(D,type="mLL",Pe=NULL,Pm=NULL,lambda=0,para=NULL,hyperpara=NULL,verbose=TRUE){				
	Sgenes = unique(colnames(D))
	n <- length(Sgenes)
	if(type %in% c("mLL","FULLmLL")) # discrete data		
		P <- sapply(Sgenes, function(x) rowSums(D[, colnames(D) == x, drop=FALSE]))
	else
		P <- D	
	cat("Estimating module network of",n,"S-genes (lambda =", lambda,")...\n\n")
	C = cov(P)
	if(all(C - diag(ncol(C)) < 1e-5)){ # if there are no correlations at all
		modeltotal <- matrix(0,ncol=n,nrow=n)	# ... there are no edges
		colnames(modeltotal) <- Sgenes
		rownames(modeltotal) <- Sgenes
	}
	else if(all(as.matrix(dist(C,method="manhattan")) < 1e-5)){ # if all genes are identical
		modeltotal <- matrix(1,ncol=n,nrow=n) # ... it has to be the fully connected graph	
		colnames(modeltotal) <- Sgenes
		rownames(modeltotal) <- Sgenes
	}	
	else{	 # otherwise something in between	
		if(is.null(Pe))
			Pe <- matrix(1/n,nrow=nrow(D),ncol=n)
		if(is.null(Pm))
			Pm <- matrix(0,ncol=n,nrow=n)		
		modeltotal <- matrix(0,ncol=n,nrow=n)
		colnames(modeltotal) <- Sgenes
		rownames(modeltotal) <- Sgenes
		modeltotal <- moduleNetwork.aux(D,modeltotal,Sgenes,Pe=Pe,Pm=Pm,lambda=lambda,para=para,hyperpara=hyperpara,type=type,verbose=verbose)
	}
	modeltotal <- transitive.closure(modeltotal, mat=TRUE,loop=TRUE)
	diag(modeltotal) <- 0
	ep <- score(list(modeltotal),D,type=type,para=para,Pe=Pe,hyperpara=hyperpara,verbose=FALSE)	
	# output
	modeltotal <- as(modeltotal,"graphNEL")	
    	res <- list(graph=modeltotal,mLL=ep$mLL[[1]],pos=ep$pos[[1]],mappos=ep$mappos[[1]],type=ep$type,para=para,hyperpara=hyperpara,lam=lambda)	
	class(res) <- "ModuleNetwork"
	if(verbose)
		cat("================ finished!!! ================\n")		
	return(res)
}

moduleNetwork.aux <- function(D,modeltotal,variables,Pe=NULL,Pm=NULL,lambda=0,para=NULL,hyperpara=NULL,type=type,verbose=TRUE){
	if(verbose){
		cat("estimating network of genes:\n")
		cat(variables,"\n")
	}
	n <- length(variables)
	if(n > 4){		
		if(type %in% c("mLL","FULLmLL")) # discrete data		
			P <- sapply(variables, function(x) rowSums(D[, colnames(D) == x, drop=FALSE]))
		else
			P <- D
		tP <- t(P)	
		C <- cov(P)
		R <- as.matrix(dist(tP,method="manhattan"))	
		if(all(C - diag(ncol(C)) < 1e-5))
			modeltotal[variables,variables] <- 0
		else if(all(R < 1e-5))
			modeltotal[variables,variables] <- 1		
		else{											
			Dis <- as.matrix(dist(tP))			
			kmax <- min(min(apply(Dis,1,function(x) length(which(x != 0)))), floor(ncol(P)/2))
			sil <- 0						
			if(verbose)
				cat(paste("--> performing clustering: can form at most ",kmax,"clusters ...\n"))	
			clusterData <- function(k, tP){
# 				pr <- kmeans(tP,k,nstart=100)
				pr <- pam(tP,k,metric="manhatten")				
# 				pr$clustering <- pr$cluster
				s <- summary(silhouette(pr,dist(tP)))$avg.width
				return(list(sil=s, clust=pr, k=k))
			}
			
			clust <- lapply(2:kmax, clusterData, tP)			
			maidx <- which.max(sapply(clust, function(x) x$sil))
			kbest <- clust[[maidx]]$k
			clust <- clust[[maidx]]$clust								
			if(verbose){
				cat(paste("---->", kbest,"clusters built\n"))				
				cat("(cluster indices:",clust$clustering,")\n")
			}
			
			for(i in 1:kbest){
				cl <- which(clust$clustering == i)							
				if(length(cl) > 1){
					Ccl = C[cl,cl,drop=FALSE]
					if(all(Ccl - diag(ncol(Ccl)) < 1e-5))
						modeltotal[variables[cl],variables[cl]] <- 0
					else if(all(R[cl,cl,drop=FALSE] < 1e-5))
						modeltotal[variables[cl],variables[cl]] <- 1					
					else{			
						vars <- colnames(P)[cl]							
						modeltotal <- moduleNetwork.aux(D[,colnames(D) %in% vars,drop=FALSE], modeltotal,variables[cl],Pe[,cl],Pm[cl,cl],lambda,para,hyperpara,type=type,verbose=verbose)
					}					
				}
			}						
			for(i in 1:kbest){
				cli <- which(clust$clustering == i)							
				if(i < kbest){
					for(j in (i+1):kbest){					
						clj <- which(clust$clustering == j)	
						if((length(cli) > 0) && (length(clj) > 0)){
							Ccl <- C[c(cli,clj),c(cli,clj),drop=FALSE]
							if(all(Ccl - diag(ncol(Ccl)) < 1e-5)){
								idx <- c(variables[cli],variables[clj])		
								modeltotal[idx,idx] <- 0
							}
							else if(all(R[c(cli,clj),c(cli,clj),drop=FALSE] < 1e-5)){
								idx <- c(variables[cli],variables[clj])		
								modeltotal[idx,idx] <- 1
							}							
							else								
								modeltotal <- connectModules(modeltotal, D, cli, clj,variables,Pe=Pe,Pm=Pm,lambda,para=para,hyperpara=hyperpara,type=type,verbose=verbose)
						}
					}
				}
			}
		}
	}
	else{		
		if(verbose)
			cat("--> estimating local model\n")
		if(n > 1){			
			models <- enumerate.models(n,name=unique(colnames(D)),verbose)
			sco <- score(models,D,type=type,para=para,hyperpara=hyperpara,Pe=Pe,Pm=Pm,lambda=lambda,verbose=verbose)$mLL
			modellocal <- models[[which.max(sco)]]			
		}
		else
			modellocal <- 0
		modeltotal[variables,variables] <- modellocal		
	}
	return(modeltotal)
}

connectModules <- function(modeltotal, D, cl1, cl2, variables, Pe=NULL, Pm=NULL,lambda=0,para=NULL,hyperpara=NULL,type="mLL",verbose=TRUE){
	if(verbose)
		cat("--> estimating connections between module",t(variables[cl1]),"and",t(variables[cl2]),"(length of module 1 = ",length(cl1),"length of module 2 = ", length(cl2),")\n")		
	idx <- c(variables[cl1],variables[cl2])		
	cl <- c(cl1,cl2)		
	modeltotalOrig <- modeltotal
	for(i in 1:length(cl1)){				
		for(j in (length(cl1)+1):length(idx)){									
# 		for(j in 1:length(cl2)){
			models <- list()
			k <- 1	
			model <- modeltotalOrig[idx,idx]
# 			Variante 1: A B unverbunden
			model <- transitive.closure(model, mat=TRUE,loop=TRUE)				
			models[[k]] <- model				
			k <- k + 1
# 			Variante 2: A -> B
			model[i,j] <- 1				
			model <- transitive.closure(model, mat=TRUE,loop=TRUE)				
			models[[k]] <- model				
			k <- k + 1
# 			Variante 3: B -> A	
			model <- modeltotalOrig[idx,idx]
			model[j,i] <- 1				
			model <- transitive.closure(model, mat=TRUE,loop=TRUE)				
			models[[k]] <- model				
			k <- k + 1
# 			Variante 4: A <-> B
			model <- modeltotalOrig[idx,idx]
			model[i,j] <- 1
			model[j,i] <- 1
			model <- transitive.closure(model, mat=TRUE,loop=TRUE)				
			models[[k]] <- model				
			k <- k + 1		
						
			bestmodel <- as(score(models,D[,colnames(D) %in% idx,drop=FALSE], type=type, para=para, hyperpara=hyperpara,Pe=Pe[,cl],Pm=Pm[cl,cl],lambda=lambda,verbose=FALSE)$graph,"matrix")
			modeltotal[idx[i],idx[j]] <- bestmodel[i,j]
			modeltotal[idx[j],idx[i]] <- bestmodel[j,i]
			
			if(verbose)
				cat(".")
		}					
	}
	if(verbose)
		cat("\n")	
	return(modeltotal)
}
