pairwise.posterior <- function(D, type="mLL", para=NULL, hyperpara=NULL,
                               Pe=NULL, Pm=NULL, verbose=TRUE) {
                               
  # Sgenes
  Sgenes <- unique(colnames(D))
  nrS <- length(Sgenes)
  nrTest <- nrS*(nrS-1)/2
  if(verbose) cat(nrS,"perturbed genes ->", nrTest, "pairwise tests\n")


  # local model prior
  if (is.null(Pm)) Pm <- rep(.25,4)

  # init output
  graph <- diag(nrS)
  dimnames(graph) <- list(Sgenes,Sgenes)
  scores <- matrix(nrow=nrTest,ncol=5)
  dimnames(scores) <- list(as.character(1:nrTest),c("..","->","<-","<->","support")) 
  ix <- 1

  # loop over edges
  for (i in 1:(nrS-1)) {
    for (j in (i+1):nrS) {

      # data
      x <- Sgenes[i]
      y <- Sgenes[j]
      D.xy <- D[ , which(colnames(D)==x | colnames(D)==y)]
      D.xy <- D.xy[rowSums(D.xy)!=0 ,,drop=FALSE]
      support <- nrow(D.xy)
      
      # four models per edge: x..y  x->y  x<-y  x<->y
      models <- enumerate.models(2,name=c(x,y),verbose=FALSE)

      # score
      ss <- score(models, D.xy, type=type, para=para,
                  hyperpara=hyperpara, Pe=Pe, verbose=FALSE)
      
      post <- exp(ss$mLL) * Pm  # model likelihood times prior
      post <- post/sum(post)    # scaled to [0,1]
      
      # winner
      winner <- models[[which.max(post)]]
      graph[i,j] <- winner[1,2]
      graph[j,i] <- winner[2,1]

      # scores
      scores[ix,] <- c(post,support)
      rownames(scores)[ix] <- paste(c(x,y),collapse=":")

      # counter :(
      ix <- ix + 1

      if(verbose) cat(".")
    }
  }
  if(verbose) cat("\n")  
  
  # estimate effect positions
  if (verbose) cat("estimating effect positions\n")
  ep <- score(list(transitive.closure(graph,mat=TRUE)), D, 
               type=type, 
               para=para,
               hyperpara=hyperpara,
               Pe=Pe, 
               verbose=FALSE)
  
  # remove loops in graph
  graph <- graph-diag(nrS)
  
  # as graphNEL
  graph <- as(graph,"graphNEL")
  
  # output
  res <- list(graph=graph,pos=ep$pos[[1]],mappos=ep$mappos[[1]],scores=scores,type=type,para=para,hyperpara=hyperpara)
  class(res) <- "pairwise"

  return(res)
}

      
