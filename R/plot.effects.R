plot.effects <- function(x,nem,border=TRUE,...){
	
int        <- unique(colnames(x))
sccg       <- SCCgraph(nem$graph,name=TRUE)
topo.order <- tsort(sccg$graph)

#----------------------------
# reorder COLUMNS            
#----------------------------

# order
ord        <- unlist(sccg$scc[topo.order])
col.order  <- unlist(lapply(ord,function(y) which(colnames(x) == y)))
D          <- x[,col.order]
Dcn        <- colnames(D)

#----------------------------
# reorder ROWS               
#----------------------------

# estimate effect positions
# old: take them directly from hierarchy:
#      if (class(nem)=="pairwise") mappos <- nem$mappos
#      if (class(nem)=="score")    mappos <- nem$mappos[[which.max(nem$mLL)]]
# new: estimate them from scc graph
# colnames(D) <- sccg$which.scc[Dcn]
# M <- as(sccg$graph,"matrix") + diag(length(sccg$scc)) 
# mappos <- score(D,models=M,type=nem$type,para=nem$para,hyperpara=nem$hyperpara,verbose=FALSE)$mappos[[1]]
# colnames(D) <- Dcn
if(length(nem$mLL) > 1)
	mappos = nem$mappos[[which.max(nem$mLL)]]
else
	mappos = nem$mappos

# order 
v <- list()
nr <- list()
for (i in topo.order){
  w = unlist(mappos[sccg$scc[[i]]], use.names=FALSE)
  
  if (length(w)==0){ 
    v[[i]]  <- NA
    nr[[i]] <- 1 
    }
  if (length(w)==1){
    v[[i]] <- w
    nr[[i]] <- 1
  }  
  if (length(w) >1){
        d       <- dist(D[w,,drop=FALSE],method="manhattan")
        v[[i]]  <- w[hclust(d)$order]
        nr[[i]] <- length(w) 
    }    
}
v <- rev(unlist(v))
nr <- rev(unlist(nr))

D2 <- matrix(0,nrow=length(v),ncol=ncol(D))
D2[which(!is.na(v)),] <- D[v[which(!is.na(v))],]
dimnames(D2) <- list(rownames(D)[v],colnames(D))

#----------------------------
# PLOT                       
#----------------------------
cs <- cumsum(nr)[-length(nr)]

par(las=2,mgp=c(5.5,1,0),mar=c(7,7,4,7),cex.lab=1.7,cex.main=2,lwd=2)
image(x   = 1:ncol(D2),
      y   = 1:nrow(D2),
      z   = t(D2),
      xlab= "",
      ylab= "",
      xaxt= "n",
      yaxt= "n",
      col = gray(seq(.95,0,length=10)),
      ...
      )
axis(1,1:ncol(D2),colnames(D2))
axis(4,1:nrow(D2),rownames(D2),tick=!is.na(v))
axis(2, c(0,cs) + nr/2 , rev(topo.order),tick=FALSE)
if (border) abline(h=cs+.5,col="red")
}
