# log marginal likelihood of model
mLL <- function(Phi,D1,D0=NULL,a=0.15,b=0.05,Pe=NULL,Pm=NULL,lambda=0,type="mLL") {       
   if (!all(diag(Phi)==1)) 
  	diag(Phi) <- 1 	        
  if(type == "mLL")
  	L  <- a^(D1 %*% (1-Phi)) * (1-a)^(D0 %*% (1-Phi)) * (1-b)^(D1 %*% Phi) * b^(D0 %*% Phi)  	    
  else if(type %in% c("CONTmLLBayes", "CONTmLLDens"))	
	L <- exp(D1%*%Phi)
  else if(type == "CONTmLL")
	L <- exp(log(D1)%*%Phi + log((1-D1))%*%(1-Phi)) 
  else if(type %in% c("CONTmLLMAP", "CONTmLLRatio")){		
	Sgenes = colnames(Phi)
	Phi = cbind(Phi, double(ncol(Phi)))
	colnames(Phi)[ncol(Phi)] = "null"
	ep = D1%*%Phi + log(Pe)	
	Theta = apply(ep,1,function(e) e ==max(e))		
	L = t((Phi%*%(Theta*1)>0)*1)*D1[,Sgenes]
	LLperGene=rowSums(L)			
	s = sum(LLperGene)	
	map = apply(Theta,1,which)	
  }
  if(!(type %in% c("CONTmLLMAP","CONTmLLRatio"))){
	if(!is.null(Pe))
		LP <- L*Pe
	else
		LP <- L	
	LLperGene = log(rowSums(LP))
	ep <- LP/(rowSums(LP))  	  			
	Theta = apply(ep,1,function(e) e ==max(e))
	s  <- sum(LLperGene)
  	map = apply(Theta,1,which)			
  }      
  if(!is.null(rownames(D1)))
	map = sapply(map, names)  
  if((lambda != 0) && !is.null(Pm))
  	s <- s - lambda*sum(abs(Phi - Pm))    
  list(mLL=s,pos=ep,mappos=map,LLperGene=LLperGene)
}

