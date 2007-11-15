nem.calcSignificance <- function(D, Phi, likelihood, N=1000, seed=1, type="mLL", Pe=NULL, para=NULL, hyperpara=NULL, delta=1, selEGenes=FALSE){
	sampleRandomNetwork = function(scaleFree=FALSE, gamma=3, maxOutDegree=3, maxInDegree=1){	
		n = ncol(Phi)
		S = diag(n) # network of S-genes	
		maxOutDegree = min(n,maxOutDegree)	
		degprob = (0:maxOutDegree)^(-gamma)
		degprob[1] = 1
		degprob = degprob/sum(degprob)	
		for(i in 1:n){ # connect network randomly
			if(scaleFree)
				outdeg = sample(0:maxOutDegree,1,prob=degprob)# power law for out-degree => scale-free network
			else
				outdeg = sample(0:maxOutDegree,1)		
			if(outdeg > 0){
				idx0 = which(S[i,] == 0)
				if(length(idx0) > 0){
					idx = which(colSums(S[,idx0]) <= maxInDegree)
					if(length(idx) > 0){			
						idx = sample(idx0[idx],min(outdeg,length(idx0[idx])),replace=TRUE)	
						S[i,idx] = 1
						S = transitive.closure(S, mat=TRUE,loop=TRUE)
					}			
				}	
			}
		}						
		diag(S) = 0			
		colnames(S) = colnames(Phi)
		rownames(S) = rownames(Phi)
		S
	}
	
	set.seed(seed)
	BFrand = double(N)
	for(i in 1:N){
		net = sampleRandomNetwork()		
		loglik <- nem(D,models=list(net),inference="search",type=type,para=para,Pe=Pe,Pm=NULL,lambda=0,delta=delta,hyperpara=hyperpara,verbose=FALSE)$mLL
		BFrand[i] = likelihood - loglik
	}		
	loglik <-  nem(D,models=list(Phi*0),inference="search",type=type,para=para,Pe=Pe,Pm=NULL,lambda=0,delta=delta,hyperpara=hyperpara,verbose=FALSE)$mLL
	BFzero = likelihood - loglik
	loglik <- nem(D,models=list((Phi>=0)*1),inference="search",type=type,para=para,Pe=Pe,Pm=NULL,lambda=0,delta=delta,hyperpara=hyperpara,verbose=FALSE)$mLL
	BFcomplete = likelihood - loglik
	p.value = length(which(BFrand <= 0)) / N
	list(p.value=p.value, BFzero=BFzero, BFcomplete=BFcomplete)
}
