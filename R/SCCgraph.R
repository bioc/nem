SCCgraph <- function(x,name=TRUE,nlength=20){

  if (!is(x, "graphNEL") & !is(x,"matrix")) stop("Wrong class of argument 'x': must be 'graphNEL' or 'matrix'")  
  if (is(x,"matrix"))  x <- as(x,"graphNEL")

  scc   <- strongComp(x)
  N     <- length(scc)
  
  # concatenate node names in same scc
  V <- as.character(1:N)
  if (name==TRUE){
    for (i in 1:N){        
        v <- paste(scc[[i]],collapse=":")
        if (nchar(v)>nlength) v <- paste(substr(v,1,nlength-3),"...",sep="")
        V[i] <- v     
    }
    names(scc) <- V
  }

  # which node is in which scc?
  which.scc <- numeric(length(nodes(x)))
  names(which.scc)<-nodes(x)
  for (i in names(scc)) which.scc[scc[[i]]] <- i


  # build scc graph
  edL <- vector("list", length = N)
  names(edL) <- V
  for (i in names(scc)){
    vv <-  scc[[i]]
    ee <- NULL
    for (j in vv) ee <- c(ee,x@edgeL[[j]]$edges)
    ee <- which.scc[ee]
    dup <- duplicated(ee)
    ee <- ee[!dup & !(ee==i)]
    edL[[i]] <- list(edges= ee)
  }
  gR <- new("graphNEL", nodes=V, edgeL=edL, edgemode="directed")

  #-----
  names(scc) <- V

  #-----
  return(list(graph=gR,scc=scc,which.scc=which.scc))
}
