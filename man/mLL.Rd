\name{mLL}
\alias{mLL}
\title{Marginal likelihood of a phenotypic hierarchy}
\description{
  computes the marginal likelihood of observed phenotypic data given a phenotypic hierarchy. 
}
\usage{
mLL(Phi, D1, D0, a, b, Pe)
}

\arguments{
  \item{Phi}{an adjacency matrix with unit main diagonal}
  \item{D1}{count matrix: phenotypes x genes. How often did we see an effect after interventions?}
  \item{D0}{count matrix: phenotypes x genes. How often did we NOT see an effect after intervention?}
  \item{a}{false positive rate: how probable is it to miss an effect?}
  \item{b}{false negative rate: how probable is it to see a spurious effect?}
  \item{Pe}{prior of effect reporter positions in the phenotypic hierarchy}
}
\details{
  It computes the marginal likelihood of a single phenotypic hierarchy.
  Usually called from within the function \code{score}.
}
\value{
  \item{mLL}{marginal likelihood of a phenotypic hierarchy}
  \item{pos}{posterior distribution of effect positions in the hierarchy}
  \item{mappos}{Maximum aposteriori estimate of effect positions}
}

\references{Markowetz F, Bloch J, Spang R, Non-transcriptional pathway features reconstructed from secondary effects of RNA interference, Bioinformatics, 2005}

\author{Florian Markowetz <URL: http://genomics.princeton.edu/~florian>}

\seealso{\code{\link{score}}, \code{\link{FULLmLL}}}
\examples{
   data("BoutrosRNAi2002")
   result <- nem(BoutrosRNAiDiscrete[,9:16],type="mLL",para=c(.15,.05))
}
\keyword{models}
