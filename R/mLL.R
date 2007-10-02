# log marginal likelihood of model
mLL <- function(Phi,D1,D0=NULL,a=0.15,b=0.05,Pe=NULL,Pm=NULL,lambda=0,type="mLL") {    
   if (!all(diag(Phi)==1)) 
  	diag(Phi) <- 1 	      
  if(type == "mLL")
  	L  <- a^(D1 %*% (1-Phi)) * (1-a)^(D0 %*% (1-Phi)) * (1-b)^(D1 %*% Phi) * b^(D0 %*% Phi)  	    
#   else if(type == "FULLCONTmLLDens")
# 	L <- fitBeta(Phi,D1,Pe)
  else if(type == "CONTmLLDens")
	L <- exp(log(D1+1e-10)%*%Phi)
  else if(type == "CONTmLL")
	L <- exp(log(D1+1e-10)%*%Phi + log((1-D1)+1e-10)%*%(1-Phi)) 
  else if(type == "CONTmLLRatio"){	
	Theta <- matrix(0,nrow=nrow(Phi),ncol=nrow(D1))			
	idx <- apply(D1%*%Phi + log(Pe),1,which.max)					
	for(r in 1:nrow(Theta))
		Theta[r,which(idx == r)] <- 1		
	L <- Phi%*%Theta%*%D1
  }
  if(type != "CONTmLLRatio"){
	if(!is.null(Pe))
		LP <- L*Pe+1e-10  
	else
		LP <- L + 1e-10
	s  <- sum(log(rowSums(LP)))
	ep <- LP/(rowSums(LP))  
  	map <- colnames(D1)[apply(ep,1,which.max)]
  	names(map) <- rownames(D1)  
  }  
  else{
	s <- sum(diag(L))
	ep <- map <- colnames(D1)[idx]
  }
  if((lambda != 0) && !is.null(Pm))
  	s <- s - lambda*sum(abs(Phi - Pm))  
  return(list(mLL=s,pos=ep,mappos=map))
}

