# log marginal likelihood of models

mLL <- function(Phi,D1,D0,a,b,Pe) {
  if (!all(diag(Phi)==1)) stop("\nnem:mLL> Model main diagonal must be 1!")
  
  # make sure model is in right order wrt to count matrices
  Phi <- Phi[colnames(D1),colnames(D1)]

  # compute the log-mll score
  L  <- a^(D1 %*% (1-Phi)) * (1-a)^(D0 %*% (1-Phi)) * (1-b)^(D1 %*% Phi) * b^(D0 %*% Phi)
  LP <- L*Pe
  s  <- sum(log(rowSums(LP)))

  # posterior of effect positions
  ep <- LP/rowSums(LP)
  colnames(ep) = colnames(D1) 
 
  # MAP estimate of effect positions
  map <- colnames(D1)[apply(ep,1,which.max)]
  names(map) <- rownames(D1)

  return(list(mLL=s,pos=ep,mappos=map))
}
