nem.greedyMAP <- function(D, Pe=NULL, Pm=NULL, lambda=0, delta=1, trans.close=TRUE, verbose=TRUE){		
	Sgenes = unique(colnames(D))
	n <- length(Sgenes)		
	cat("Alternating optimization for",n,"S-genes (lambda =", lambda,")...\n\n")
	Theta = apply(D,1, function(e) (e == max(e)) & (e>0))*1			
	best = NULL
	converged = FALSE
	i = 1		
	while(!converged){
# 		if(trans.close & (n < 9)){
# 			cvec = as.vector(t(Theta%*%D))									
# 			T = diag(n)			
# 			combis = combn(which(T==0),3)
# 			bvec = rep(1,ncol(combis)+n)
# 			Amat = matrix(0,nrow=ncol(combis),ncol=length(cvec))
# 			for(c in 1:ncol(combis)){
# 				Amat[c,combis[,c]] = c(1,1,-1)
# 			}					
# 			ones = which(T==1)			
# 			for(c in 1:length(ones)){				
# 				tmp = double(n^2)
# 				tmp[ones[c]] = 1
# 				Amat = rbind(Amat, tmp)
# 			} 
# 			
# 			Phi0 = lp("max", cvec, Amat, c(rep("<=",ncol(combis)),rep("=",n)), bvec, all.bin=TRUE)		
# 			Phi = matrix(Phi0$solution, byrow=FALSE, ncol=n)			
# 			dimnames(Phi) = list(Sgenes, Sgenes)	
# 			print(Phi)
# 			readline()					
# 		}
# 		else
			Phi = apply(Theta%*%D, 1, function(e) (e > 0))*1
		sco <- score(list(Phi),D,type="CONTmLLMAP",Pe=Pe,Pm=Pm,lambda=lambda,delta=delta,verbose=verbose)
		if(verbose){
			cat("iteration ", i, ": likelihood = ", sco$mLL, "\n")
		}
		if(i > 1)		
			converged = all(Phi==Phiold)
		if(converged)
			break;	
		Theta = apply(sco$pos[[1]],1,function(e) (e == max(e)))*1				
		Theta = Theta[1:(nrow(Theta)-1),]				
		Phiold = Phi		
		i = i + 1
	}			
	Phi = as(sco$graph, "matrix")	
	if(trans.close)	
		Phi = closest.transitive.greedy(Phi, verbose)
	ep <- score(list(Phi),D,type="CONTmLLMAP",Pe=Pe,Pm=Pm,lambda=lambda,delta=delta,verbose=verbose)
    	res <- list(graph=ep$graph,mLL=ep$mLL[[1]],pos=ep$pos[[1]],mappos=ep$mappos[[1]],type=ep$type[[1]],para=ep$para,hyperpara=ep$hyperpara,lam=lambda,selected=ep$selected,delta=delta,LLperGene=ep$LLperGene[[1]])	# output: data likelihood under given model!	
	class(res) <- "nem.greedyMAP"	
	return(res)
}
