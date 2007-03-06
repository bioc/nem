\name{FULLmLL}
\alias{FULLmLL}
\title{Full marginal likelihood of a phenotypic hierarchy}
\description{
  The function the full marginal likelihood of a phenotypic hierarchy.
  The full marginal likelihood equals the marginal likelihood \code{mLL} 
  averaged over the error probabilities $\alpha$ and $\beta$.
  
}
\usage{
FULLmLL(Phi, D1, D0, a0, b0, a1, b1, Pe, Pm=NULL, lambda=0)
}

\arguments{
  \item{Phi}{an adjacency matrix with unit main diagonal}
  \item{D1}{count matrix: phenotypes x genes. How often did we see an effect after interventions?}
  \item{D0}{count matrix: phenotypes x genes. How often did we NOT see an effect after intervention?}
  \item{a0, b0, a1, b1}{Hyperparameters}
  \item{Pe}{prior of effect positions in the hierarchy. A matrix of size phenotypes x genes, where each row contains positive numbers summing to 1.}
  \item{Pm}{prior on model graph (n x n matrix) with entries 0 <= priorPhi[i,j] <= 1 describing the probability of an edge between gene i and gene j.}
  \item{lambda}{regularization parameter to incorporate prior assumptions.}
}
\details{
  Additionally to the marginal likelihood introduced in Markowetz et al (2005),
  we can average over the error probabilities $\alpha$ and $\beta$ assuming Beta priors.
  The parameters of the two Beta priors are hyperparameters of the full marginal likelihood score.
  The four hyperparameters fall into two categories: 
  \code{a1} and \code{b0} are weights for observing the predicted state, while \code{a0} and \code{b1} 
  are weights for observing errors.
  We suggest setting \code{a1=b0} and \code{a0=b1}. 
  The ratio between the two values should correspond to our assessment of the noise level.
  See the example section for an application. 
  The function \code{FULLmLL} is usually called from within function \code{score}.
}
\value{
  \item{mLL}{full marginal likelihood of a model}
  \item{pos}{posterior distribution of effect positions in the hierarchy}
  \item{mappos}{maximum aposteriori estimate of effect positions}

}

\references{Markowetz F, Probabilistic Models for Gene Silencing
  Data. PhD thesis, Free University Berlin, 2006.}
\author{Florian Markowetz <URL: http://genomics.princeton.edu/~florian>}

\seealso{\code{\link{score}}, \code{\link{mLL}}}
\examples{   
   data("BoutrosRNAi2002")
   res <- nem(BoutrosRNAiDiscrete[,9:16],type="FULLmLL",hyperpara=c(1,9,9,1))
}
\keyword{models}
