\name{nem.greedyMAP}
\alias{nem.greedyMAP}
\alias{print.nem.greedyMAP}

\title{Infers a phenotypic hierarchy using an alternating MAP optimization}
\description{
   Starting with an initial estimate of the linking of E-genes to S-genes from the data, this method performs an alternating MAP optimization of the S-genes graph and the linking graph until convergence. As a final step the function \code{closest.transitive.greedy} can be invoked to find a transitively closed graph most similar to the original one.
}
\usage{
nem.greedyMAP(D,Pe=NULL,Pm=NULL,lambda=0,delta=1, trans.close=TRUE, verbose=TRUE)

\method{print}{nem.greedyMAP}(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{  
\item{D}{data matrix. Columns correspond to the nodes in the silencing scheme. Rows are phenotypes.}  
  \item{Pe}{prior position of effect reporters. Default: uniform over nodes in hierarchy}  
  \item{Pm}{prior on model graph (n x n matrix) with entries 0 <= priorPhi[i,j] <= 1 describing the probability of an edge between gene i and gene j.}
  \item{lambda}{regularization parameter to incorporate prior assumptions.}  
  \item{delta}{regularization parameter for automated E-gene subset selection}
  \item{trans.close}{find a similar transitively closed graph}
  \item{verbose}{do you want to see progress statements printed or not? Default: TRUE}  


  \item{x}{nem object}
  \item{...}{other arguments to pass}
}
\value{
  nem object
}
\author{Holger Froehlich}

\examples{
   data("BoutrosRNAi2002") 
   res <- nem(BoutrosRNAiLods, inference="nem.greedyMAP", delta=0)
   
   # plot graph
   plot(res,what="graph")
   
   # plot posterior over effect positions
   plot(res,what="pos")
   
   # estimate of effect positions
   res$mappos
   
}

\seealso{\code{\link{nem}}, \code{\link{closest.transitive.greedy}}}

\keyword{models}
