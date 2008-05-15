\name{pairwise.posterior}
\alias{pairwise.posterior}
\alias{print.pairwise}

\title{Infers a phenotypic hierarchy edge by edge}
\description{
  Function \code{pairwise.posterior} estimates the hierarchy edge by edge. In each step only a pair of nodes
  is involved and no exhaustive enumeration of model space is needed as in function \code{score}.
}
\usage{
pairwise.posterior(D, type = "mLL", para = NULL, hyperpara = NULL, 
    Pe = NULL, Pmlocal = NULL, Pm = NULL, lambda = 0, delta=1, verbose = TRUE)

\method{print}{pairwise}(x,...) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{D}{data matrix. Columns correspond to the nodes in the silencing scheme. Rows are phenotypes.}  
\item{type}{see \code{nem}}    
  \item{para}{vector with parameters a and b for "mLL", if count matrices are used}
  \item{hyperpara}{vector with hyperparameters a0, b0, a1, b1 for "FULLmLL"}
  \item{Pe}{prior position of effect reporters. Default: uniform over nodes in hierarchy}
  \item{Pmlocal}{local model prior for the four models tested at each node: a vector of length 4 with positive entries summing to one}
  \item{Pm}{prior on model graph (n x n matrix) with entries 0 <= priorPhi[i,j] <= 1 describing the probability of an edge between gene i and gene j.}
  \item{lambda}{regularization parameter to incorporate prior assumptions.}  
  \item{delta}{regularization parameter for automated E-gene subset selection (CONTmLLRatio only)}  
  \item{verbose}{do you want to see progress statements printed or not? Default: TRUE} 

  \item{x}{nem object}
  \item{...}{other arguments to pass}
}
\details{
  \code{pairwise.posterior} is a fast(er) heuristic alternative to exhaustive search 
  by the function \code{score}. For each pair (\code{A},\code{B}) of perturbed genes 
  it chooses between four possible models: \code{A..B} (unconnected), \code{A->B} (superset), 
  \code{A<-B} (subset), or \code{A<->B} (indistinguishable). 
  The result is the graph built from the maximum aposteriori models for each edge.
  
 \code{print.pairwise} gives an overview over the 'pairwise' object.
}

\value{
nem object
}

\author{Florian Markowetz <URL: http://genomics.princeton.edu/~florian>}


\seealso{\code{\link{score}}, \code{\link{nem}}}
\examples{
   data("BoutrosRNAi2002") 
   res <- nem(BoutrosRNAiDiscrete[,9:16],para=c(.13,.05),inference="pairwise")
   
   # plot graph
   plot(res,what="graph")
   
   # plot posterior over effect positions
   plot(res,what="pos")
   
   # estimate of effect positions
   res$mappos
   
}
\keyword{graphs}
