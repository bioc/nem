triples.posterior <- function(D,type="mLL",para=NULL, hyperpara=NULL,Pe=NULL,Pm=NULL,triples.thrsh=.5,verbose=TRUE){

  # Sgenes
  Sgenes <- unique(colnames(D))
  nrS    <- length(Sgenes)
  nrTest <- choose(nrS,3) 
  if(verbose) cat(nrS,"perturbed genes ->", nrTest, "triples to check\n")

  # local model prior
  if (is.null(Pm)) Pm <- rep(1/29,29)
 
 
  ## 
  ## 1. learn tripel models
  ##
 
  ##=== create all (ordered) triples
  triples = subsets(nrS,3)

  ##=== 29 candidate models for each triple
  cnd.models   = enumerate.models(3,verbose=FALSE)

  ##=== store maximum model in a list
  mll.models   = list()
  if (verbose) prev <- progressBar()
  for(i in 1:nrow(triples)){

        D.tmp = D[,colnames(D)%in%Sgenes[triples[i,]]]
        D.tmp <- D.tmp[rowSums(D.tmp)!=0,]
        colnames(D.tmp)[colnames(D.tmp) == Sgenes[triples[i,1]]] = "a"
        colnames(D.tmp)[colnames(D.tmp) == Sgenes[triples[i,2]]] = "b"
        colnames(D.tmp)[colnames(D.tmp) == Sgenes[triples[i,3]]] = "c"
        tmp.mdl = score(cnd.models,D.tmp,para=para,hyperpara=hyperpara, Pe=Pe, verbose=FALSE)
 
        ##== do prior stuff
        post = tmp.mdl$mLL + log(Pm) ##=== prior needs to be in enumerate.models(3) order!
        winner <- cnd.models[[which.max(post)]]
        mll.models[[i]] = list()
        mll.models[[i]]$graph = winner
        mll.models[[i]]$posterior = post
        dimnames(mll.models[[i]]$graph) = list(Sgenes[triples[i,]],Sgenes[triples[i,]]) #=== rename node to inds
        if (verbose) prev <- progressBar(i/nrTest,prev)
  }

  ##
  ## 2. break down to pairwise-data by model averaging
  ##

  A = matrix(NA, ncol=nrS, nrow=nrS)
  dimnames(A) = list(Sgenes,Sgenes)

  diag(A) = 1
  for(i in 1:(nrS-1)){  
  for(j in (i+1):nrS){
        ##=== find contributing triples
        inds = apply(triples,1,function(x) sum( x %in% c(i,j)) == 2)
        inds = which(inds)
        ##=== count edges and set A_ij accordingly
        tmp = list()
        for(k in 1:length(inds)){ ##=== for each contributing triple
                tmp[[k]] = mll.models[[inds[k]]]$graph #=== get adjacency matrix
        }
        ##=== mean number of edges from i to j (including doubles...)
        A[i,j] = mean(unlist(lapply(tmp,function(x) x[Sgenes[i],Sgenes[j]])))
        ##=== mean number of edges from j to i (including doubles...) 
        A[j,i] = mean(unlist(lapply(tmp,function(x) x[Sgenes[j],Sgenes[i]])))
  }}
  B = (A>=triples.thrsh)*1-diag(nrow(A))
  graph <- as(B,"graphNEL")



  ##
  ## 3. estimate effect positions
  ##
  if (verbose) cat("Estimating effect positions in combined graph\n")
  ep <- score(list(transitive.closure(graph,mat=TRUE)), D, 
               type=type, 
               para=para,
               hyperpara=hyperpara,
               Pe=Pe, 
               verbose=FALSE)

  ##
  ## 4. output
  ##
  res <- list(graph=graph,avg=A,pos=ep$pos[[1]],mappos=ep$mappos[[1]],type=type,para=para,hyperpara=hyperpara)
  class(res) <- "triples"
  
  return(res)
}