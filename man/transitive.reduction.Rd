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
  \code{transitive.reduction} uses a modification of the classical algorithm from the Sedgewick book for computing transitive closures.
}
\value{
returns a graph object with shortcuts removed
}
\references{
R. Sedgewick, Algorithms, Pearson, 2002.
}

\author{Holger Froehlich}

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
