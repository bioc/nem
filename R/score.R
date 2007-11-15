score <- function(models, D, type="mLL", para=NULL, hyperpara=NULL, Pe=NULL, Pm=NULL, lambda=0, delta=1, verbose=TRUE, graphClass="graphNEL") {

  #if single model as input
  if (class(models)=="matrix") models <- list(models)    

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
  if (!all(sapply(models,fkt,s=Sgenes))) stop("\nnem:score> models must have same names as data")
  nrS <- length(Sgenes)     
  # make probability/density matrices D0 and D1  
  # nrow=#E-genes and ncol=#S-genes        
  if(type %in% c("CONTmLL","CONTmLLRatio","CONTmLLDens")){   		
  	# D1[i,j] = probability/density of EFFECT at E_i when S_j was silenced  	  	  		
	D1 <- D
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
    # if no prior is supplied:
  # assume uniform prior over E-gene positions      
  if (is.null(Pe)){ 	
	Pe <- matrix(1/nrS,nrow=nrow(D1),ncol=nrS)
	colnames(Pe) <- Sgenes  		
  }          
  if(type == "CONTmLLRatio"){		
  	Pe = cbind(Pe, double(nrow(D1)))  			
  	Pe[,ncol(Pe)] = delta/nrS
	Pe = Pe/rowSums(Pe)	
  }
  if(is.null(Pm)) lambda <- 0 
   
    
  if (type=="FULLmLL"){ # FULL log marginal likelihood of all models
    if (verbose==TRUE) cat("Computing FULL marginal likelihood for",length(models),"models\n")
    a0 <- hyperpara[1]
    b0 <- hyperpara[2]
    a1 <- hyperpara[3]
    b1 <- hyperpara[4]    
    results <- sapply(models,FULLmLL,D1,D0,a0,b0,a1,b1,Pe=Pe,Pm=Pm,lambda=lambda)        
  }
  else{   # log marginal likelihood of all models	
	if (verbose==TRUE) cat("Computing marginal likelihood for",length(models),"models\n")      	
	a <- para[1]
	b <- para[2]               
	results <- sapply(models,mLL,D1,D0,a,b,Pe=Pe,Pm=Pm,lambda=lambda,type=type)     	   			
  }     
  s       <- unlist(results["mLL",])
  ep      <- results["pos",]
  map     <- results["mappos",]   
  if(!is.null(Pm)){  	
  	log_pD_cond_Phi <- s  	  	  		
  	if(is.null(lambda) || (lambda == 0)){  			
		if(verbose) cat("--> Using Bayesian model averaging to incorporate prior knowledge\n")
  		lpPhi <- sapply(models, PhiDistr, Pm, a=1, b=0.5)		
  		ppost <- exp((log_pD_cond_Phi + lpPhi)/nrS)	
  		s <- log_pD_cond_Phi + lpPhi  		
  	}
  	else
  		ppost <- exp((s + log(lambda*0.5))/nrS)
  }  
  else
	ppost = exp(s/nrS)
  ppost <- ppost/sum(ppost)  	# posterior model probability

  
  if(verbose){
	if(length(ppost) > 1){
		sort_ppost <- sort(ppost,decreasing=TRUE)
		cat("(probabilities of best and second best model for ",Sgenes, ":", sort_ppost[1],",",sort_ppost[2],")\n")
	}
  }  
   # winning model       
  winner <- models[[which.max(s)]]  
  diag(winner) <- 0  
  if(graphClass == "graphNEL"){
  	gR <- new("graphAM",adjMat=winner,edgemode="directed")  
  	gR <- as(gR,"graphNEL")    
  }
  else
	gR <- winner
  PHI <- matrix(0,ncol=nrS,nrow=nrS)
  dimnames(PHI) <- list(Sgenes,Sgenes)  
  for(i in 1:length(models))
	PHI <- PHI + models[[i]]*ppost[i]           
  selected = which(unlist(results["mappos",which.max(s)]) != "NA")	    
  # output  
  res <- list(graph=gR, mLL=s, ppost=ppost, avg=PHI, pos=ep, mappos=map, type=type, para=para, hyperpara=hyperpara, lam=lambda, selected=selected)
  class(res) <- "score"   
  return(res)  
}
