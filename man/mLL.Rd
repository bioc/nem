\name{mLL}
\alias{mLL}
\title{Marginal likelihood of a phenotypic hierarchy}
\description{
  computes the marginal likelihood of observed phenotypic data given a phenotypic hierarchy. 
}
\usage{
mLL(Phi,D1,D0=NULL,a=0.15,b=0.05,Pe=NULL,Pm=NULL,lambda=0,type="mLL")
}

\arguments{
  \item{Phi}{an adjacency matrix with unit main diagonal}
  \item{D1}{(i) count matrix for discrete data: phenotypes x genes. How often did we see an effect after interventions? (ii) matrix describing the PROBABILITIES of an effect (iii) matrix discribing the log-LIKELIHOOD of an effect (e.g. log-density matrix, log-odds matrix)}
  \item{D0}{count matrix: phenotypes x genes. How often did we NOT see an effect after intervention? Not used for continious data}
  \item{a}{false positive rate: how probable is it to miss an effect? (for count matrix)}
  \item{b}{false negative rate: how probable is it to see a spurious effect? (for count matrix)}
  \item{Pe}{prior of effect reporter positions in the phenotypic hierarchy}
  \item{Pm}{prior on model graph (n x n matrix) with entries 0 <= priorPhi[i,j] <= 1 describing the probability of an edge between gene i and gene j.}
  \item{lambda}{regularization parameter to incorporate prior assumptions.}
  \item{type}{see \code{nem}}
}
\details{
  It computes the marginal likelihood of a single phenotypic hierarchy.
  Usually called from within the function \code{score}.
}
\value{
  \item{mLL}{marginal likelihood of a phenotypic hierarchy}
  \item{pos}{posterior distribution of effect positions in the hierarchy}
  \item{mappos}{Maximum aposteriori estimate of effect positions}
  \item{LLperGene}{likelihood per E-gene}
}

\references{Markowetz F, Bloch J, Spang R, Non-transcriptional pathway features reconstructed from secondary effects of RNA interference, Bioinformatics, 2005}

\author{Holger Froehlich <URL: http://www.dkfz.de/mga2/people/froehlich>, Florian Markowetz <URL: http://genomics.princeton.edu/~florian>}

\seealso{\code{\link{nem}}, \code{\link{score}}, \code{\link{FULLmLL}}}
\examples{
   data("BoutrosRNAi2002")
   result <- nem(BoutrosRNAiDiscrete[,9:16],type="mLL",para=c(.15,.05))
}
\keyword{models}
