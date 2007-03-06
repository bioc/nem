transitive.reduction <- function(g){
	if(class(g) == "matrix"){
		# modified algorithm from Sedgewick book: just remove transitive edges instead of inserting them
		g = g - diag(diag(g))
		for(y in 1:nrow(g)){
			for(x in 1:nrow(g)){
				if(g[x,y] != 0){
					for(j in 1:nrow(g)){
						if(g[y,j] != 0)
							g[x,j] = 0
					}
				}
			}
		}
	}
	else{
		with.children <- sapply(edgeL(g),function(x) length(x$edges)>0)
		# loop over nodes with children
		for(i in nodes(g)[with.children]){
			visited <- rep(FALSE,length(nodes(g)))
			names(visited) <- nodes(g)
			#loop over children of 'i' which have children of their own
			for (j in names(which(with.children[adj(g,i)[[1]]])) ){
			# loop over grandchildren of 'i' which have not been visited yet
				for (k in names(which(!visited[adj(g,j)[[1]]]))){
				# if grandchild can also be reached directly -> remove it
					if (k %in% adj(g,i)[[1]]){
					g <- removeEdge(i,k,g)
					visited[k] <- TRUE
				}
			}
			}
		}
	}
	return(g)
}