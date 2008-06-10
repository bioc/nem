plotnem = function(D, G, x, SCC, main=NULL, zlim=NULL){
	if(length(x$mLL) > 1){
		winner <- which.max(x$mLL)
		mappos <- x$mappos[[winner]]		
	}
	else{
		mappos <- x$mappos		
	}						
	nf = split.screen(rbind(c(0,1,0.4,1),c(0,0.8,0,0.4), c(0.9,1,0,0.4)))
	erase.screen(2)
	screen(2)	
	mynodes = AgNode(G)
	nodenames = sapply(mynodes, name)		
	xy = getNodeXY(G)
	left = xy$x <= max(xy$x)*0.5
	right = xy$x > max(xy$x)*0.5	
	nodenames = c(nodenames[left][order(xy$x[left] + xy$y[left])], nodenames[right][order(xy$x[right] + xy$y[right],decreasing=TRUE)])		
	if(!is.null(zlim))
		ord = plotEffects(D, x, legend=FALSE, order=nodenames, orderSCC=SCC, zlim=zlim)		
	else
		ord = plotEffects(D, x, legend=FALSE, order=nodenames, orderSCC=SCC)		
	erase.screen(3)		
	screen(3)
	nrcolors = 200; half = 1+nrcolors/2 # nrcolors must be an even number
	colpal = c(brewer.pal(9,"Blues")[9:1],brewer.pal(9,"OrRd")[1:9])
	allcolors = colorRampPalette(colpal)(nrcolors)
	if(is.null(zlim)){
		r = c(quantile(D[D<0],0.95), quantile(D[D>0],0.95))
		rangeall = c(-max(abs(r),na.rm=TRUE), max(abs(r),na.rm=TRUE))	
	}
	else
		rangeall = zlim
	color.legend(0.4,0.1,1,1,signif(seq(rangeall[1],rangeall[2],length.out=5),digits=1),rect.col=allcolors[length(allcolors):1],gradient="y", cex=0.75)

	erase.screen(1)
	screen(1)
	par(mar=c(0,0,0,0))
	plot(G, main=main)				
	ma = 0.82*getX(upRight(boundBox(G)))		
	xrescale = function(x,a2=0,b2=ma){
		alpha = (a2-b2)/(1-length(ord))
		beta = a2 - alpha*1		
		alpha*x + beta
	}
	for(i in 1:length(mynodes)){
		xy = getNodeXY(mynodes[[i]])		
		nodenames = strsplit(name(mynodes[[i]]),":")[[1]]
		to = xrescale(match(unique(unlist(mappos[nodenames])),ord))
		for(j in 1:length(to)){			
			lines(xy.coords(c(xy$x, to[j]),c(xy$y-0.5*getNodeHeight(mynodes[[i]]), 0)), col="grey", lty=1, lwd=1, type="b", pch=19, cex=0.3)	
		}		
	}	
	screen(1)
	plot(G, main=main)	
	close.screen(c(1,2,3),all.screens=TRUE)					
}