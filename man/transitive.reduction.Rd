\name{transitive.reduction}
\alias{transitive.reduction}

\title{Computes the transitive reduction of a graph}
\description{
  \code{transitive.reduction} removes direct edges, which can be explained by another path in the graph.
}
\usage{
transitive.reduction(g)
}

\arguments{
  \item{g}{graphNEL object}
}
\details{
  \code{transitive.reduction} implements an interative algorithm in spirit of Wagner (2001).
}
\value{
returns a graph object with shortcuts removed
}
\references{Wagner A, How to reconstruct a large genetic network from $n$ gene perturbations in fewer than $n^2$ easy steps. Bioinformatics 17(12): 1183-1197 (2001)}
\author{Florian Markowetz <URL: http://genomics.princeton.edu/~florian>}


\seealso{\code{\link{transitive.closure}}}
\examples{
   V <- LETTERS[1:3]
   edL <- list(A=list(edges=c("B","C")),B=list(edges="C"),C=list(edges=NULL))
   gc <- new("graphNEL",nodes=V,edgeL=edL,edgemode="directed")
   g <- transitive.reduction(gc)
    
   par(mfrow=c(1,2))
   plot(gc,main="shortcut A->C")
   plot(g,main="shortcut removed")
}
\keyword{graphs}
