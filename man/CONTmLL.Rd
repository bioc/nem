\name{CONTmLL}
\alias{CONTmLL}
\title{Marginal likelihood of a phenotypic hierarchy with continuous data}
\description{
  computes the marginal likelihood of observed phenotypic data given a phenotypic hierarchy.
}
\usage{
CONTmLL(Phi, prob.inf, Pe)
}

\arguments{
  \item{Phi}{an adjacency matrix with unit main diagonal}
  \item{prob.inf}{matrix of influence probabilities}
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

\references{ }

\author{Matthias Maneck}

\seealso{\code{\link{score}}, \code{\link{mLL}}, \code{\link{FULLmLL}}, \code{\link{nem.cont.preprocess}}}
\examples{
   data("BoutrosRNAi2002")
   pre <- nem.cont.preprocess(BoutrosRNAiExpression,neg.control=1:4,pos.control=5:8)
   result <- nem(pre$prob.influenced, type="CONTmLL")
}
\keyword{models}



