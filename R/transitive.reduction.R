transitive.reduction <-
function(g){
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
return(g)
}
