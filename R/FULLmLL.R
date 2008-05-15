# FULL log marginal likelihood of a models

FULLmLL <- function(Phi,D1,D0,a0,b0,a1,b1,Pe,Pm=NULL,lambda=0) {
  if (!all(diag(Phi)==1)) diag(Phi) = 1
  
  # make sure model is in right order wrt to count matrices
  Phi <- Phi[colnames(D1),colnames(D1)]

  # compute score
  n01 <- D1 %*% (1-Phi)
  n00 <- D0 %*% (1-Phi)
  n11 <- D1 %*% Phi
  n10 <- D0 %*% Phi
  s0  <- gamma(a0+b0)*gamma(n10+a0)*gamma(n00+b0)/gamma(a0)/gamma(b0)/gamma(n10+n00+a0+b0)
  s1  <- gamma(a1+b1)*gamma(n11+a1)*gamma(n01+b1)/gamma(a1)/gamma(b1)/gamma(n11+n01+a1+b1)  
  SP  <- s0*s1*Pe
  LLperGene = log(rowSums(SP))
  s   <- sum(LLperGene)
  
  if((lambda != 0) && !is.null(Pm))
  	s <- s - lambda*sum(abs(Phi - Pm))

  # posterior effect positions
  ep  <- SP/rowSums(SP)

  # MAP estimate of effect positions
  Theta = apply(ep,1,function(e) e ==max(e))
  map = apply(Theta,1,which)	

  return(list(mLL=s,pos=ep,mappos=map,LLperGene=LLperGene))
}
