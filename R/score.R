score <- function(models,D,type="mLL", para=NULL, hyperpara=NULL, Pe=NULL,verbose=TRUE) {

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

  # if no prior is supplied:
  # assume uniform prior over E-gene positions
  if (is.null(Pe)) Pe <- matrix(1/nrS,nrow=nrow(D),ncol=nrS)
  
  # make count matrices D0 and D1
  # nrow=#E-genes and ncol=#S-genes
  # D0[i,j] = how often NO EFFECT at E_i when S_j was silenced
  # D1[i,j] = how often    EFFECT at E_i when S_j was silenced
  D0  <- matrix(NA,ncol=nrS,nrow=nrow(D),dimnames=list(rownames(D),Sgenes))
  D1  <- D0
  for (i in 1:nrS) {
    Di     <- D[,colnames(D) == Sgenes[i],drop=FALSE]
    D0[,i] <- rowSums(Di==0)
    D1[,i] <- rowSums(Di==1)    
  }  

  # log marginal likelihood of all models
  if (type=="mLL"){
    if (verbose==TRUE) cat("Computing marginal likelihood for",length(models),"models\n")
    a <- para[1]
    b <- para[2]
    results <- sapply(models,mLL,D1,D0,a,b,Pe)
    s       <- unlist(results["mLL",])
    ep      <- results["pos",]
    map     <- results["mappos",]
  }

  # FULL log marginal likelihood of all models
  if (type=="FULLmLL"){
    if (verbose==TRUE) cat("Computing FULL marginal likelihood for",length(models),"models\n")
    a0 <- hyperpara[1]
    b0 <- hyperpara[2]
    a1 <- hyperpara[3]
    b1 <- hyperpara[4]
    results <- sapply(models,FULLmLL,D1,D0,a0,b0,a1,b1,Pe)
    s       <- unlist(results["mLL",])
    ep      <- results["pos",]
    map     <- results["mappos",]

  }

  # winning model
  winner <- as(models[[which.max(s)]],"graphNEL")

  # output
  res <- list(graph=winner,mLL=s,pos=ep, mappos=map, type=type, para=para, hyperpara=hyperpara)
  class(res) <- "score"
  return(res)  
}
