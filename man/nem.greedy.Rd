\name{nem.greedy}
\alias{nem.greedy}
\alias{print.nem.greedy}

\title{Infers a phenotypic hierarchy using a greedy search strategy}
\description{
   Starting from an initial graph (default: no edges), this strategy successively adds those edges, which most inrease the likelihood of the data under the model.
}
\usage{
nem.greedy(D,initial=NULL,type="mLL",Pe=NULL,Pm=NULL,lambda=0,delta=1,para=NULL,hyperpara=NULL,verbose=TRUE)

\method{print}{nem.greedy}(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{  
\item{D}{data matrix. Columns correspond to the nodes in the silencing scheme. Rows are phenotypes.}  
\item{initial}{initial model to start greedy hillclimbing from (default: no edges)} 
\item{type}{(1.) marginal likelihood "mLL" (only for cout matrix D), or (2.) full marginal likelihood "FULLmLL" integrated over a and b and depending on hyperparameters a0, a1, b0, b1 (only for count matrix D), or (3.) "CONTmLL" marginal likelihood for probability matrices, or (4.) "CONTmLLDens" marginal likelihood for probability log-density matrices, or (5.) "CONTmLLRatio" for log-odds ratio matrices}  
  \item{Pe}{prior position of effect reporters. Default: uniform over nodes in hierarchy}  
  \item{Pm}{prior on model graph (n x n matrix) with entries 0 <= priorPhi[i,j] <= 1 describing the probability of an edge between gene i and gene j.}
  \item{lambda}{regularization parameter to incorporate prior assumptions.}  
  \item{delta}{regularization parameter for automated E-gene subset selection (CONTmLLRatio only)}
  \item{para}{vector with parameters a and b for "mLL", if count matrices are used}
  \item{hyperpara}{vector with hyperparameters a0, b0, a1, b1 for "FULLmLL"}  
  \item{verbose}{do you want to see progress statements printed or not? Default: TRUE}  


  \item{x}{nem object}
  \item{...}{other arguments to pass}
}
\value{
  \item{graph}{the inferred directed graph (graphNEL object)}   
  \item{mLL}{marginal likelihood of final model} 
  \item{pos}{posterior over effect positions}
  \item{mappos}{MAP estimate of effect positions} 
  \item{type}{as used in function call}
  \item{para}{as used in function call}
  \item{hyperpara}{as used in function call}
  \item{lambda}{as in function call}
  \item{selected}{selected E-gene subset}
}
\author{Holger Froehlich}

\seealso{\code{\link{nem}}}
\examples{
   # Drosophila RNAi and Microarray Data from Boutros et al, 2002
   data("BoutrosRNAi2002")
   D <- BoutrosRNAiDiscrete[,9:16]
   nem.greedy(D, para=c(.13,.05))
}

\keyword{models}
