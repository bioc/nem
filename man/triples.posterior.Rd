\name{triples.posterior}
\alias{triples.posterior}
\alias{print.triples}

\title{Infers a phenotypic hierarchy from triples}
\description{
  Function \code{triples.posterior} estimates the hierarchy triple-wise. In each step only a triple of nodes
  is involved and no exhaustive enumeration of model space is needed as in function \code{score}.
}
\usage{
triples.posterior(D, type="mLL",para=NULL, hyperpara=NULL,Pe=NULL,Pmlocal=NULL,Pm=NULL,lambda=0,delta=1, triples.thrsh=.5,verbose=TRUE)

\method{print}{triples}(x,...) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{D}{data matrix. Columns correspond to the nodes in the silencing scheme. Rows are phenotypes.}   
\item{type}{(1.) marginal likelihood "mLL" (only for cout matrix D), or (2.) full marginal likelihood "FULLmLL" integrated over a and b and depending on hyperparameters a0, a1, b0, b1 (only for count matrix D), or (3.) "CONTmLL" marginal likelihood for probability matrices, or (4.) "CONTmLLDens" marginal likelihood for probability density matrices, or (5.) "CONTmLLRatio" for log-odds ratio matrices}  
  \item{para}{vector with parameters a and b for "mLL", if count matrices are used}
  \item{hyperpara}{vector with hyperparameters a0, b0, a1, b1 for "FULLmLL"}
  \item{Pe}{prior position of effect reporters. Default: uniform over nodes in hierarchy}
  \item{Pmlocal}{local model prior for the four models tested at each node: a vector of length 4 with positive entries summing to one}
  \item{triples.thrsh}{threshold used when combining tripel models for each edge. Default: only edges appearing in more than half of triples are included in the final graph.}
  \item{Pm}{prior on model graph (n x n matrix) with entries 0 <= priorPhi[i,j] <= 1 describing the probability of an edge between gene i and gene j.}
  \item{lambda}{regularization parameter to incorporate prior assumptions.}  
  \item{delta}{regularization parameter for automated E-gene subset selection (CONTmLLRatio only)}
  \item{verbose}{do you want to see progress statements printed or not? Default: TRUE}  

  \item{x}{nem object}
  \item{...}{other arguments to pass}
}

\details{
  \code{triples.posterior} is an alternative to exhaustive search 
  by the function \code{score} and more accurate than \code{pairwise.posterior}. 
  For each triple of perturbed genes 
  it chooses between the 29 possible models. It then uses model averaging to combine the triple-models into a final graph.
  
  
  \code{print.triples} gives an overview over the 'triples' object.
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

\references{
  Markowetz F, Kostka D, Troyanskaya OG, Spang R:
  Nested effects models for high-dimensional phenotyping screens.
  Bioinformatics. 2007; 23(13):i305-12.
}

\author{Florian Markowetz <URL: http://genomics.princeton.edu/~florian>}


\seealso{\code{\link{score}}, \code{\link{nem}}}
\examples{
   data("BoutrosRNAi2002") 
   res <- triples.posterior(BoutrosRNAiDiscrete[,9:16],para=c(.13,.05))
   
   # plot graph
   plot(res,what="graph")
   
   # plot posterior over effect positions
   plot(res,what="pos")
   
   # estimate of effect positions
   res$mappos
   
}
\keyword{graphs}
