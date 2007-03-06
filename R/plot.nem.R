plot.nem <- function(x, what="graph", remove.singletons=FALSE, PDF=FALSE, filename="nemplot.pdf", thresh=0, transitiveReduction=FALSE, ...) {
	
    if (!(what%in%c("graph","mLL","pos"))) stop("\nnem> invalid plotting type: plot either 'graph', 'mLL', or 'pos'")

    if (what=="graph"){
    	if(class(x) != "matrix")
        	M <- as(x$graph,"matrix")
        else
        	M <- x
        if (all(diag(M)==1)) M <- M-diag(ncol(M))
        if (remove.singletons){
            take  <- colSums(M) != 0 | rowSums(M) != 0
            M <- M[take,take]
        }        
        if(thresh != 0){ # discretize matrix based on threshold  			
		nodes <- colnames(M)
		nodenames <- vector("character",length(nodes)^2)
		probs <- vector("double",length(nodes)^2)
		k <- 1
		for(i in 1:length(nodes)){
			for(j in 1:length(nodes)){
				probs[k] <- as.character(M[i,j])
				nodenames[k] <- paste(nodes[i],"~",nodes[j],sep="")
				k <- k + 1
			}  	
		}  	
		names(probs) <- nodenames
		probs <- probs[probs > thresh]
		edgeattr <- list(label=probs)  	  	
		M[M <= thresh] <- 0  	  	
  	}  
        if(transitiveReduction==TRUE)
  		M <- transitive.reduction(M)  
	gR <- new("graphAM",adjMat=M,edgemode="directed")	
	if (PDF) pdf(file=filename)   
        par(cex.main=2)        
        if(thresh == 0)
		plot(x=gR, y="dot",...)
	else
		plot(x=gR, y="dot",edgeAttrs=edgeattr,...)                
        if (PDF) dev.off()
    }

    if(what=="mLL"){
       par(cex=1.3)
      ss <- sort(unique(x$mLL),decreasing=TRUE)[1:min(30,length(x$mLL))]
      plot(x=1:length(ss), y=ss, pch=19, main="Score distribution",
        xlab=paste(length(ss),"top ranked models"),
        ylab="Marginal log-likelihood", 
        ylim=c(ss[length(ss)]-10,ss[1]+10)
        )
      points(1,max(unique(x$mLL)),pch=21,cex=1.7,lwd=2)

    }
    
    if(what=="pos"){    
	if(length(x$mLL) > 1){
        	winner <- which.max(x$mLL)
		pos <- x$pos[[winner]]
		effects <- rownames(x$pos[[winner]])
	}
	else{
		pos <- x$pos
		effects <- rownames(x$pos)
	}
        par(las=2,mgp=c(5.5,1,0),mar=c(6.7,7,4,1),cex.lab=1.3,cex.main=1.7)        
        image(x=1:4,
            y=1:nrow(pos),
            z = t(pos),
            main = "Posterior effect positions",
            xlab="Perturbations",
            xaxt="n",
            ylab="Effect reporters",
            yaxt="n",
            col=gray(seq(.95,0,length=10))
        )
        abline(v=(1:3)+.5)
        axis(1,1:4,colnames(x$graph))        
        axis(2,1:length(effects),effects)
    }

}
  
