sampleRndNetwork = function(Sgenes, scaleFree=TRUE, gamma=3, maxOutDegree=3, maxInDegree=2, trans.close=TRUE){		
	n = length(Sgenes)
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
				}			
			}	
		}
	}						
	if(trans.close)
		S = transitive.closure(S, mat=TRUE,loop=FALSE)					
	diag(S) = 0			
	colnames(S) = Sgenes
	rownames(S) = Sgenes
	S
}

sampleData = function(Phi, m, prob=NULL, uninformative=0, type="binary", alpha=sample(seq(0.1,0.9,by=0.1),ncol(Phi),replace=TRUE), beta=sample(5:50,ncol(Phi),replace=TRUE), lambda=matrix(sample(seq(0.01,0.49,by=0.01),ncol(Phi)*2,replace=TRUE),ncol=2), meansH1=rep(0.5, ncol(Phi)), meansH0=rep(-0.5, ncol(Phi)), sdsH1=sample(seq(0.1,1,by=0.1),ncol(Phi),replace=TRUE), sdsH0=sample(seq(0.1,1,by=0.1),ncol(Phi),replace=TRUE)){
	Sgenes = colnames(Phi)
	n = length(Sgenes)
	epos = sample(1:n,m,replace=TRUE,prob=prob)
	Theta = matrix(0, nrow=length(Sgenes), ncol=m)
	for(i in 1:ncol(Theta))
		Theta[epos[i],i] = 1
	M = t(Phi%*%Theta)
	if(type == "binary"){
		D = M
		if(uninformative > 0)
			D = rbind(D, matrix(0,ncol=n, nrow=uninformative))			
	}
	else if(type %in% c("density")){
		lambda = cbind(1-rowSums(lambda), lambda)				
		palt = sapply(1:n, function(i) bum.ralt(m, c(alpha[i], beta[i]), lambda[i,]))		
		p0 = matrix(runif((m+uninformative)*n), ncol=n)
		P = M*palt + (1-M)*p0[1:m,]
		if(uninformative > 0)
			P = rbind(P, p0[(m+1):nrow(p0),])
		D = sapply(1:n, function(i) bum.dalt(P[,i], c(alpha[i], beta[i]), lambda[i,]))
		D = log(D)
	}
	else if(type %in% c("lodds")){
		palt = sapply(1:n, function(i) rnorm(m, mean=meansH1[i], sd=sdsH1[i]))
		p0 = sapply(1:n, function(i) rnorm(m+uninformative, mean=meansH0[i], sd=sdsH0[i]))		
		D = M*palt + (1-M)*p0[1:m,]
		if(uninformative > 0)
			D = rbind(D, p0[(m+1):nrow(p0),])		
	}	
	else
		stop(paste("unknown type", type, "\n"))	
	colnames(D) = Sgenes
	D
}
