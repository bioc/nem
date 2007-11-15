moduleNetwork <- function(D,type="mLL",Pe=NULL,Pm=NULL,lambda=0,delta=1,para=NULL,hyperpara=NULL,verbose=TRUE){				
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
				modeltotal <- moduleNetwork.aux(D[,varidx,drop=FALSE],modeltotal,variables,Pe=Pe[,variables],Pm=Pm[variables,variables],lambda=lambda,delta=delta,para=para,hyperpara=hyperpara,type=type,verbose=verbose)
			}	
		}		
	}
	if(length(idx) > 1)
		res <- connectModules(D, modeltotal, idx, Pe=Pe,Pm=Pm,lambda=lambda,delta=delta,para=para,hyperpara=hyperpara,type=type,verbose=verbose)
	else{
		modeltotal <- transitive.closure(modeltotal, mat=TRUE,loop=TRUE)	
		diag(modeltotal) <- 0	
		ep <- score(list(modeltotal),D,type=type,para=para,Pe=Pe,Pm=NULL,lambda=0,delta=delta,hyperpara=hyperpara,verbose=FALSE)		
	# output
		modeltotal <- as(modeltotal,"graphNEL")	
		res <- list(graph=modeltotal,mLL=ep$mLL[[1]],pos=ep$pos[[1]],mappos=ep$mappos[[1]],type=ep$type,para=para,hyperpara=hyperpara,lam=lambda,selected=ep$selected)	# output: data likelihood under given model!
		class(res) <- "ModuleNetwork"
		if(verbose){
			cat("log-likelihood of model = ",res$mLL,"\n")
		}
	}
	return(res)
}

moduleNetwork.aux <- function(D,modeltotal, variables,Pe=NULL,Pm=NULL,lambda=0,delta=1,para=NULL,hyperpara=NULL,type="mLL",verbose=TRUE){
	if(verbose){
		cat("estimating network of genes:\n")
		cat(variables,"\n")
	}
	n <- length(variables)
	if(n > 1){			
		models <- enumerate.models(n,name=unique(colnames(D)),verbose=verbose)
		sco <- score(models,D,type=type,para=para,hyperpara=hyperpara,Pe=Pe,Pm=Pm,lambda=lambda,delta=delta,verbose=verbose,graphClass="matrix")
		modellocal <- sco$graph
	}
	else
		modellocal <- 0
	modeltotal[variables,variables] <- modellocal
	modeltotal
}

connectModules <- function(D, Phi, modules, type="mLL",Pe=NULL,Pm=NULL,lambda=0,delta=1,para=NULL,hyperpara=NULL,verbose=TRUE){
	Sgenes = unique(colnames(D))
	n <- length(Sgenes)	
	if(verbose){	
		dev.off()
		cat("Connecting modules using constraint greedy hillclimbing ...\n\n")	
	}
	sco0 <- score(list(Phi),D,type=type,para=para,hyperpara=hyperpara,Pe=Pe,Pm=Pm,lambda=lambda,delta=delta,verbose=verbose,graphClass="matrix")$mLL			
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
			sconew <- score(models,D,type=type,para=para,hyperpara=hyperpara,Pe=Pe,Pm=Pm,lambda=lambda,delta=delta,verbose=verbose,graphClass="matrix")			
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
	ep <- score(list(Phi),D,type=type,para=para,Pe=Pe,Pm=NULL,lambda=0,delta=delta,hyperpara=hyperpara,verbose=FALSE,graphClass="matrix")
	diag(Phi) <- 0
	Phi <- as(Phi,"graphNEL")	
    	res <- list(graph=Phi,mLL=ep$mLL[[1]],pos=ep$pos[[1]],mappos=ep$mappos[[1]],type=ep$type,para=para,hyperpara=hyperpara,lam=lambda,selected=ep$selected)	# output: data likelihood under given model!
	class(res) <- "ModuleNetwork"
	if(verbose)
		cat("log-likelihood of model = ",res$mLL,"\n")
	return(res)
}

# moduleNetwork.auxOld <- function(D,modeltotal,idx, variables,Pe=NULL,Pm=NULL,lambda=0,para=NULL,hyperpara=NULL,type="mLL",selEGenes=FALSE,verbose=TRUE){
# 	if(verbose){
# 		cat("estimating network of genes:\n")
# 		cat(variables,"\n")
# 	}
# 	n <- length(variables)
# 	if(n > 4){				
# 		P <- sapply(variables, function(x) rowSums(D[, colnames(D) == x, drop=FALSE]))			
# 		tP <- t(P)	
# 		C <- abs(cov(P))
# 		diag(C) = 0
# 		R <- as.matrix(dist(tP,method="manhattan"))	
# 		if(all(C < 1e-5)){
# 			modeltotal[variables,variables] <- 0
# 		}
# 		else if(all(R < 1e-5)){
# 			modeltotal[variables,variables] <- 1		
# 		}
# 		else{						
# 			Dis0 <- dist(tP,method="manhattan")								
# 			Dis <- as.matrix(Dis0)			
# 			kmax <- min(min(apply(Dis,1,function(x) length(which(x != 0))+1)), floor(ncol(P)/2))
# 			sil <- 0						
# 			if(verbose)
# 				cat(paste("--> performing clustering: can form at most ",kmax,"clusters ...\n"))
# 			hc <- hclust(Dis0,method="ward")
# 			clusterData <- function(k, tP){
# 				pr <- list(clustering=cutree(hc,k=k))
# 				s <- summary(silhouette(pr,Dis0))$avg.width
# 				return(list(sil=s, clust=pr, k=k))
# 			}
# 			
# 			clust <- lapply(2:kmax, clusterData, tP)			
# 			maidx <- which.max(sapply(clust, function(x) x$sil))
# 			kbest <- clust[[maidx]]$k
# 			clust <- clust[[maidx]]$clust								
# 			if(verbose){
# 				cat(paste("---->", kbest,"clusters built\n"))				
# 				cat("(cluster indices:",clust$clustering,")\n")
# 			}
# 			
# 			for(i in 1:kbest){
# 				cl <- which(clust$clustering == i)							
# 				if(length(cl) > 1){					
# 					if(all(C[cl,cl,drop=FALSE] < 1e-5)){
# 						modeltotal[variables[cl],variables[cl]] <- 0
# 					}
# 					else if(all(R[cl,cl,drop=FALSE] < 1e-5)){
# 						modeltotal[variables[cl],variables[cl]] <- 1			
# 					}		
# 					else{										
# 						vars <- colnames(P)[cl]							
# 						varidx <- sapply(vars,function(x) which(colnames(D) %in% x))		
# 						result <- moduleNetwork.aux(D[,varidx,drop=FALSE], modeltotal,idx, variables[cl],Pe[,cl],Pm[cl,cl],lambda,para,hyperpara,type=type,selEGenes=selEGenes,verbose=verbose)
# 						modeltotal <- result$model
# 						idx <- result$idx
# 					}					
# 				}
# 				else
# 					idx = c(idx, list(variables[cl]))
# 			}						
# 		}
# 	}
# 	else{		
# 		if(verbose)
# 			cat("--> estimating local model\n")
# 		if(n > 1){			
# 			models <- enumerate.models(n,name=unique(colnames(D)),verbose)
# 			sco <- score(models,D,type=type,para=para,hyperpara=hyperpara,Pe=Pe,Pm=Pm,lambda=lambda,selEGenes=selEGenes,verbose=verbose)
# 			modellocal <- as(sco$graph,"matrix") 			
# 		}
# 		else
# 			modellocal <- 0
# 		modeltotal[variables,variables] <- modellocal
# 		idx = c(idx,list(variables))		
# 	}
# 	return(list(model=modeltotal,idx=idx))
# }
# 
# connectModulesOld <- function(modeltotal, D, cl1, cl2, variables, Pe=NULL, Pm=NULL,lambda=0,para=NULL,hyperpara=NULL,type="mLL",selEGenes=FALSE,verbose=TRUE){
# 	if(verbose)
# 		cat("--> estimating connections between module",t(variables[cl1]),"and",t(variables[cl2]),"(length of module 1 = ",length(cl1),"length of module 2 = ", length(cl2),")\n")				
# # 	ntriples = choose(length(cl11),2)*length(cl22)	
# # 	if(ntriples > 200){					
# 		if(verbose)
# 			cat(" using pairwise inference\n")
# 		idx <- c(variables[cl1],variables[cl2])				
# 		varidx <- sapply(idx,function(x) which(colnames(D) %in% x))
# 		cl <- c(cl1,cl2)
# 		modeltotalOrig <- modeltotal		
# 		for(i in 1:length(cl1)){				
# 			for(j in (length(cl1)+1):length(idx)){	
# 				models <- list()
# 				k <- 1	
# 				model <- modeltotalOrig[idx,idx]
# 	# 			Variante 1: A B unverbunden
# 				model <- transitive.closure(model, mat=TRUE,loop=TRUE)				
# 				models[[k]] <- model				
# 				k <- k + 1
# 	# 			Variante 2: A -> B
# 				model[i,j] <- 1				
# 				model <- transitive.closure(model, mat=TRUE,loop=TRUE)				
# 				models[[k]] <- model				
# 				k <- k + 1
# 	# 			Variante 3: B -> A	
# 				model <- modeltotalOrig[idx,idx]
# 				model[j,i] <- 1				
# 				model <- transitive.closure(model, mat=TRUE,loop=TRUE)				
# 				models[[k]] <- model				
# 				k <- k + 1
# 	# 			Variante 4: A <-> B
# 				model <- modeltotalOrig[idx,idx]
# 				model[i,j] <- 1
# 				model[j,i] <- 1
# 				model <- transitive.closure(model, mat=TRUE,loop=TRUE)				
# 				models[[k]] <- model				
# 				k <- k + 1	
# 		
# 				if(verbose)
# 					cat("Connection ",idx[i], "-", idx[j], ":)\n")	
# 					
# 				sco <- score(models,D[,varidx,drop=FALSE], type=type, para=para, hyperpara=hyperpara,Pe=Pe[,cl],Pm=Pm[cl,cl],lambda=lambda,selEGenes=selEGenes,verbose=verbose)
# 				
# 				bestmodel <- as(sco$graph,"matrix")
# 				modeltotal[idx[i],idx[j]] <- bestmodel[i,j]
# 				modeltotal[idx[j],idx[i]] <- bestmodel[j,i]						
# 			}		
# 		}			
# # 		idx <- c(variables[cl1],variables[cl2])	
# # 		cl <- c(cl1,cl2)		
# # 		for(i in 1:length(cl1)){				
# # 			for(j in (length(cl1)+1):length(idx)){	
# # 				models <- enumerate.models(2,name=c(idx[i],idx[j]),verbose)
# # 				if(verbose)
# # 					cat("Connection ",idx[i], "-", idx[j], ":)\n")	
# # 				sel <- which(colnames(D) %in% c(idx[i],idx[j]))            
# #             			D.xy <- D[, sel]
# # 				sel <- c(cl[i], cl[j])
# # 				Pesel <- Pe[, sel, drop=FALSE]             			
# #             			Pmsel <- Pm[sel, sel, drop=FALSE]
# # 				sco <- score(models, D.xy, type = type, para = para,
# # 				hyperpara = hyperpara, Pe = Pesel, Pm=Pmsel, lambda=lambda, selEGenes=selEGenes, verbose = FALSE)
# # 				bestmodel <- as(sco$graph,"matrix")
# # 				modeltotal[idx[i],idx[j]] <- bestmodel[1,2]
# # 				modeltotal[idx[j],idx[i]] <- bestmodel[2,1]						
# # 			}
# # 		}	
# # 	}
# # 	else{
# # 		if(length(cl1) >= length(cl2)){
# # 			cl11 <- cl1
# # 			cl22 <- cl2
# # 		}
# # 		else{
# # 			cl11 <- cl2
# # 			cl22 <- cl1
# # 		}	
# # # 		if(verbose)
# # # 			cat(" using triple inference (", ntriples, "triples )\n")
# # 		idx <- c(variables[cl11],variables[cl22])
# # 		cl <- c(cl11,cl22)
# # 		triples <- subsets(length(cl11),2)
# # 		cnd.models <- enumerate.models(3,verbose=FALSE)		
# # 		mll.models <- list()
# # 		k <- 1
# # 		for(i in 1:nrow(triples)){
# # 			cl1sel <-  idx[triples[i,]]
# # 			for(j in (length(cl11)+1):length(idx)){
# # 				cl2sel <- idx[j]				
# # 				sel <- c(cl1sel,cl2sel)
# # 				Pesel <- Pe[,sel]				
# # 				Pmsel <- Pm[sel,sel]					
# # 				varidx <- sapply(sel,function(x) which(colnames(D) %in% x))
# # 				D.tmp <- D[,varidx,drop=FALSE]					 	
# #         			# restrict model space by including what we already know => just 16 models        	
# #         			model <- modeltotal[sel,sel]      					
# #         			cnd.models2 <- unique(lapply(cnd.models,
# #         				function(x){ 
# #         					xx <- model
# #         					xx[cl1sel[1],cl2sel] <- x["a","c"]
# #         					xx[cl2sel,cl1sel[1]] <- x["c","a"]
# #         					xx[cl1sel[2],cl2sel] <- x["b","c"]
# #         					xx[cl2sel,cl1sel[2]] <- x["c","b"]
# # 						return(xx)
# #         				}))           				
# #         			sco <- score(cnd.models2, D.tmp, type=type, para=para,hyperpara=hyperpara, Pe=Pesel, Pm=Pmsel, lambda=lambda,selEGenes=selEGenes,verbose=FALSE)   			    
# # 				bestmodel <- as(sco$graph,"matrix")							
# # 				mll.models[[k]] <- bestmodel								
# # 				k <- k + 1				
# # 				if(verbose)
# # 					cat(".")
# # 			}
# # 		}										
# #   		for(i in 1:length(cl11)){  
# #   			for(j in (length(cl11)+1):length(idx)){
# #   				sel <- c(idx[i],idx[j])
# # 				contrvars <- sapply(mll.models, function(x) colnames(x))
# # 				inds <- which(apply(contrvars, 2, function(x) all(sel %in% x)))				
# # 				##=== mean number of edges from i to j (including doubles...)
# # 				tmp <- list()
# # 				for(k in 1:length(inds))
# # 					tmp[[k]] <- mll.models[[inds[k]]]					
# # 				modeltotal[idx[i],idx[j]] = (mean(sapply(tmp,function(x) x[idx[i],idx[j]])) >= 0.5)*1		
# #   			}
# #   		}  	  		
# # 	}
# 	if(verbose)
# 		cat("\n")		
# 	return(modeltotal)
# }
