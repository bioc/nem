nem.greedy <- function(D,initial=NULL,type="mLL",Pe=NULL,Pm=NULL,lambda=0,delta=1,para=NULL,hyperpara=NULL,verbose=TRUE){
	Sgenes = unique(colnames(D))
	n <- length(Sgenes)		
	cat("Greedy hillclimber for",n,"S-genes (lambda =", lambda,")...\n\n")
	if(is.null(initial))
		Phi <- matrix(0,nrow=n,ncol=n)		
	else
		Phi = initial		
	diag(Phi) <- 1	
	dimnames(Phi) <- list(Sgenes,Sgenes)	
	if(verbose & !is.null(initial)){
		cat("initial network:\n")
		print(Phi)
	}
	sco0 <- score(list(Phi),D,type=type,para=para,hyperpara=hyperpara,Pe=Pe,Pm=Pm,lambda=lambda,delta=delta,verbose=verbose)$mLL	
	finished <- FALSE
	while(!finished){
# 		propose new edges		
		idx = which(Phi == 0)
		if(length(idx) > 0){
			models <- list()
			for(i in 1:length(idx)){ # test all possible new edges
				Phinew = Phi
				Phinew[idx[i]] = 1
				Phinew = transitive.closure(Phinew, mat=TRUE,loop=TRUE)	
				models[[i]] <- Phinew
			}
			models <- unique(models)
			sconew <- score(models,D,type=type,para=para,hyperpara=hyperpara,Pe=Pe,Pm=Pm,lambda=lambda,delta=delta,verbose=verbose)		
			if(max(sconew$mLL) > sco0){
				sco0 <- max(sconew$mLL)
				Phi <- as(sconew$graph,"matrix")			
			}
			else # otherwise no improving edge could be inserted
				finished <- TRUE
		}else
			finished <- TRUE	
	}
	ep <- score(list(Phi),D,type=type,para=para,Pe=Pe,Pm=NULL,lambda=0,delta=delta,hyperpara=hyperpara,verbose=FALSE)
	diag(Phi) <- 0
	Phi <- as(Phi,"graphNEL")	
    	res <- list(graph=Phi,mLL=ep$mLL[[1]],pos=ep$pos[[1]],mappos=ep$mappos[[1]],type=ep$type,para=para,hyperpara=hyperpara,lam=lambda,selected=ep$selected)	# output: data likelihood under given model!	
	class(res) <- "nem.greedy"
	if(verbose)
		cat("log-likelihood of model = ",res$mLL,"\n")
	return(res)
}

