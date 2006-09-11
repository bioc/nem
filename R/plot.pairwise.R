plot.pairwise <- function(x, what="graph", remove.singletons=FALSE, PDF=FALSE, filename="nemplot.pdf", ...) {

    if (!(what%in%c("graph","pos"))) stop("\nnem> plot either 'graph' or 'pos'")

    if (what=="graph"){
        M <- as(x$graph,"matrix")
        if (all(diag(M)==1)) M <- M-diag(ncol(M))
        if (remove.singletons){
            take  <- colSums(M) != 0 | rowSums(M) != 0
            graph <- graph[take,take]
        }
        gR <- as(M,"graphNEL")
        if (PDF) pdf(file=filename)   
        par(cex.main=2)
        plot(x=gR, y="dot", ...)
        if (PDF) dev.off()
    }

    if(what=="pos"){    
        par(las=2,mgp=c(5.5,1,0),mar=c(6.7,7,4,1),cex.lab=1.7,cex.main=2)
        pos <- x$pos
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
        effects <- rownames(x$pos)
        axis(2,1:length(effects),effects)
    }


}
  
