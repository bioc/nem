\name{triples.posterior}
\alias{triples.posterior}
\alias{plot.triples}
\alias{print.triples}

\title{Infers a phenotypic hierarchy from triples}
\description{
  Function \code{pairwise.posterior} estimates the hierarchy edge by edge. In each step only a pair of nodes
  is involved and no exhaustive enumeration of model space is needed as in function \code{score}.
}
\usage{
triples.posterior(D, type="mLL", para=NULL, hyperpara=NULL,Pe=NULL, Pm=NULL,triples.thrsh=.5, verbose=TRUE)

#S3 methods for class 'triples'
plot.triples(x, what="graph", remove.singletons=FALSE, PDF=FALSE, filename="nemplot.pdf", ...)
print.triples(x,...) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{D}{data matrix. Columns correspond to the nodes in the silencing scheme. Rows are phenotypes.}
  \item{type}{(1.) marginal likelihood "mLL" depending on paramters a and b, or (2.) full marginal likelihood "FULLmLL" integrated over a and b and depending on hyperparameters a0, a1, b0, b1 }
  \item{para}{vector with parameters a and b }
  \item{hyperpara}{vector with hyperparameters a0, b0, a1, b1}
  \item{Pe}{prior position of effect reporters. Default: uniform over nodes in hierarchy}
  \item{Pm}{local model prior for the four models tested at each node: a vector of length 4 with positive entries summing to one}
  \item{triples.thrsh}{threshold used when combining tripel models for each edge. Default: only edges appearing in more than half of triples are included in the final graph.}
  \item{verbose}{do you want to see progress statements printed or not? Default: TRUE}
  \item{x}{an object of class 'pairwise'}
  \item{what}{type of plot: 'graph' or 'pos'. Default: 'graph'}
  \item{remove.singletons}{remove single nodes which are not connected to any other node when plotting? Default: FALSE}
  \item{PDF}{output as pdf file? Default: FALSE}
  \item{filename}{name of the pdf if any. Default: "nemplot.pdf"}
  \item{...}{additional arguments for plotting}
}

\details{
  \code{triples.posterior} is an alternative to exhaustive search 
  by the function \code{score} and more accurate than \code{pairwise.posterior}. 
  For each triple of perturbed genes 
  it chooses between the 29 possible models. It then uses model averaging to combine the triple-models into a final graph.
  
  
  \code{plot.triples} plots the inferred phenotypic hierarchy as a directed graph, and \code{print.triples} gives an overview over the 'triples' object.
}

\value{
  \item{graph}{the inferred directed graph (graphNEL object)}
  \item{avg}{matrix of edge frequencies in triple models}
  \item{pos}{posterior over effect positions}
  \item{mappos}{MAP estimate of effect positions} 
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
