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
    Pe = NULL, Pmlocal = NULL, Pm = NULL, lambda = 0, selEGenes=FALSE, verbose = TRUE)

#S3 methods for class 'pairwise'
print.pairwise(x,...) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{D}{data matrix. Columns correspond to the nodes in the silencing scheme. Rows are phenotypes.}  
\item{type}{(1.) marginal likelihood "mLL" (only for cout matrix D), or (2.) full marginal likelihood "FULLmLL" integrated over a and b and depending on hyperparameters a0, a1, b0, b1 (only for count matrix D), or (3.) "CONTmLL" marginal likelihood for probability matrices, or (4.) "CONTmLLDens" marginal likelihood for probability density matrices, or (5.) "CONTmLLRatio" for log-odds ratio matrices}    
  \item{para}{vector with parameters a and b for "mLL", if count matrices are used}
  \item{hyperpara}{vector with hyperparameters a0, b0, a1, b1 for "FULLmLL"}
  \item{Pe}{prior position of effect reporters. Default: uniform over nodes in hierarchy}
  \item{Pmlocal}{local model prior for the four models tested at each node: a vector of length 4 with positive entries summing to one}
  \item{Pm}{prior on model graph (n x n matrix) with entries 0 <= priorPhi[i,j] <= 1 describing the probability of an edge between gene i and gene j.}
  \item{lambda}{regularization parameter to incorporate prior assumptions.}  
  \item{selEGenes}{optimize selection of E-genes for each model}
  \item{verbose}{do you want to see progress statements printed or not? Default: TRUE} 

  \item{x}{nem object}
  \item{...}{other arguments to pass}
}
\details{
  \code{pairwise.posterior} is a fast(er) heuristic alternative to exhaustive search 
  by the function \code{score}. For each pair (\code{A},\code{B}) of perturbed genes 
  it chooses between four possible models: \code{A..B} (unconnected), \code{A->B} (superset), 
  \code{A<-B} (subset), or \code{A<->B} (undistinguishable). 
  The result is the graph built from the maximum aposteriori models for each edge.
  
 \code{print.pairwise} gives an overview over the 'pairwise' object.
}

\value{
  \item{graph}{the inferred directed graph (graphNEL object)}
  \item{pos}{posterior over effect positions}
  \item{mappos}{MAP estimate of effect positions} 
  \item{scores}{a matrix with the posterior probabilities for each local model as rows}
  \item{type}{as used in function call}
  \item{para}{as used in function call}
  \item{hyperpara}{as used in function call}
  \item{lambda}{as in function call}
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
