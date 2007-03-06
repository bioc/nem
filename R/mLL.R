# log marginal likelihood of model
mLL <- function(Phi,D1,D0=NULL,a=0.05,b=0.15,Pe=NULL,Pm=NULL,lambda=0,type="mLL") {    
   Phi <- matrix(Phi,ncol=ncol(D1),nrow=ncol(D1))      
   if (!all(diag(Phi)==1)) 
  	diag(Phi) <- 1 	  	
  if(type == "mLL")
  	L  <- a^(D1 %*% (1-Phi)) * (1-a)^(D0 %*% (1-Phi)) * (1-b)^(D1 %*% Phi) * b^(D0 %*% Phi)  	  
  else if(type == "CONTmLLDens")
  	L <- exp(log(D1)%*%Phi) 
  else if(type == "CONTmLL")
	L <- exp(log(D1)%*%Phi + log(1-D1)%*%(1-Phi)) 
    
  LP <- L*Pe+1e-10  
  s  <- sum(log(rowSums(LP)))   
  if((lambda != 0) && !is.null(Pm))
  	s <- s - lambda*sum(abs(Phi - Pm))  
  ep <- LP/(rowSums(LP))  
  map <- colnames(D1)[apply(ep,1,which.max)]
  names(map) <- rownames(D1)  
  return(list(mLL=s,pos=ep,mappos=map))
}

