\name{pairwise.posterior}
\alias{pairwise.posterior}
\alias{plot.pairwise}
\alias{print.pairwise}

\title{Infers a phenotypic hierarchy edge by edge}
\description{
  Function \code{pairwise.posterior} estimates the hierarchy edge by edge. In each step only a pair of nodes
  is involved and no exhaustive enumeration of model space is needed as in function \code{score}.
}
\usage{
pairwise.posterior(D, type="mLL", para=NULL, hyperpara=NULL,Pe=NULL, Pm=NULL, verbose=TRUE)

#S3 methods for class 'pairwise'
plot.pairwise(x, what="graph", remove.singletons=FALSE, PDF=FALSE, filename="nemplot.pdf", ...)
print.pairwise(x,...) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{D}{data matrix. Columns correspond to the nodes in the silencing scheme. Rows are phenotypes.}
  \item{type}{(1.) marginal likelihood "mLL" depending on paramters a and b, or (2.) full marginal likelihood "FULLmLL" integrated over a and b and depending on hyperparameters a0, a1, b0, b1 }
  \item{para}{vector with parameters a and b }
  \item{hyperpara}{vector with hyperparameters a0, b0, a1, b1}
  \item{Pe}{prior position of effect reporters. Default: uniform over nodes in hierarchy}
  \item{Pm}{local model prior for the four models tested at each node: a vector of length 4 with positive entries summing to one}
  \item{verbose}{do you want to see progress statements printed or not? Default: TRUE}
  \item{x}{an object of class 'pairwise'}
  \item{what}{type of plot: 'graph' or 'pos'. Default: 'graph'}
  \item{remove.singletons}{remove single nodes which are not connected to any other node when plotting? Default: FALSE}
  \item{PDF}{output as pdf file? Default: FALSE}
  \item{filename}{name of the pdf if any. Default: "nemplot.pdf"}
  \item{...}{additional arguments for plotting}

}
\details{
  \code{pairwise.posterior} is a fast(er) heuristic alternative to exhaustive search 
  by the function \code{score}. For each pair (\code{A},\code{B}) of perturbed genes 
  it chooses between four possible models: \code{A..B} (unconnected), \code{A->B} (superset), 
  \code{A<-B} (subset), or \code{A<->B} (undistinguishable). 
  The result is the graph built from the maximum aposteriori models for each edge.
  
  \code{plot.pairwise} plots the inferred phenotypic hierarchy as a directed graph, and \code{print.pairwise} gives an overview over the 'pairwise' object.
}

\value{
  \item{graph}{the inferred directed graph (graphNEL object)}
  \item{pos}{posterior over effect positions}
  \item{mappos}{MAP estimate of effect positions} 
  \item{scores}{a matrix with the posterior probabilities for each local model as rows}
  \item{type}{as used in function call}
  \item{para}{as used in function call}
  \item{hyperpara}{as used in function call}
}

\author{Florian Markowetz <URL: http://genomics.princeton.edu/~florian>}


\seealso{\code{\link{score}}, \code{\link{nem}}}
\examples{
   data("BoutrosRNAi2002") 
   res <- pairwise.posterior(BoutrosRNAiDiscrete[,9:16],para=c(.13,.05))
   
   # plot graph
   plot(res,what="graph")
   
   # plot posterior over effect positions
   plot(res,what="pos")
   
   # estimate of effect positions
   res$mappos
   
}
\keyword{graphs}
