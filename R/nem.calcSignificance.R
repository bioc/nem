nem.calcSignificance <- function(D, x, N=1000, seed=1, Pe=NULL, Pm=NULL, selEGenes=(length(x$selected) != nrow(D))){	
	modify.rand = function(Phi0){			
		Phired = transitive.reduction(Phi0)
		actions = c("insertion","deletion")				
		if(all(Phi0 - diag(ncol(Phi0)) == 0))
			actions = "insertion"	
		else if(!any(Phi0) == 0)
			actions = setdiff(actions,"insertion")
		if(length(actions) > 1)
			act = sample(actions,1)	
		else
			act = actions[1]
		if(act == "insertion"){		
			Phinew = Phi0
			idx = sample(which(Phi0 == 0),1)
			Phinew[idx] = 1
			Phinew = transitive.closure(Phinew, mat=TRUE,loop=TRUE)					
		}
		else if(act == "deletion"){
			Phinew = Phired
			diag(Phired)=0
			idx = sample(which(Phired == 1),1)
			Phinew[idx] = 0
			Phinew = transitive.closure(Phinew, mat=TRUE,loop=TRUE)			
		}	
		return(Phinew)
	}	

	cat("Testing significance of NEM model\n")
	likelihood = x$mLL
	Phi = as(x$graph, "matrix")
	Phi[Phi != 0] = 1	
	Sgenes = unique(colnames(D))
	set.seed(seed)
			
	BFperm = double(N)	
	BFrand = double(N)	
	BFmod = double(N)
	for(i in 1:N){
		if(i%%10 == 0)
			cat(".")
		net = sampleRndNetwork(Sgenes)		
		loglik <- nem(D,models=list(net),inference="search",type=x$type,para=x$para,Pe=Pe,Pm=Pm,lambda=x$lam,delta=x$delta,hyperpara=x$hyperpara, selEGenes=selEGenes, verbose=FALSE)$mLL
		BFrand[i] = likelihood - loglik	

		rnd = sample(1:ncol(Phi), replace=FALSE)
		net = Phi[rnd,rnd]		
		dimnames(net) = list(Sgenes, Sgenes)		
		loglik <- nem(D,models=list(net),inference="search",type=x$type,para=x$para,Pe=Pe,Pm=Pm,lambda=x$lam,delta=x$delta,hyperpara=x$hyperpara, selEGenes=selEGenes, verbose=FALSE)$mLL
		BFperm[i] = likelihood - loglik	

		net = modify.rand(Phi)		
		loglik <- nem(D,models=list(net),inference="search",type=x$type,para=x$para,Pe=Pe,Pm=Pm,lambda=x$lam,delta=x$delta,hyperpara=x$hyperpara,selEGenes=selEGenes, verbose=FALSE)$mLL
		BFmod[i] = likelihood - loglik
	}		
	p.value.rnd = length(which(BFrand <= 0)) / N
	p.value.perm = length(which(BFperm <= 0)) / N
	p.value.mod = length(which(BFmod <= 0)) / N
	cat("\ndone.\n")
	
	list(p.value.rnd=p.value.rnd, p.value.perm=p.value.perm, p.value.mod=p.value.mod)
}
