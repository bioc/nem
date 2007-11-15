\name{nem.calcSignificance}
\alias{nem.calcSignificance}

\title{Statistical significance of network hypotheses}
\description{   
 Assess statistical significance of a network hypothesis compared to a random, an empty and a fully connected network. 
}
\usage{
	nem.calcSignificance(D, Phi, likelihood, N=1000, seed=1, type="mLL", Pe=NULL, para=NULL, hyperpara=NULL, delta=1, selEGenes=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{D}{data matrix with experiments in the columns (binary or continious)}  
  \item{Phi}{network hypothesis (adjacency matrix)}
  \item{likelihood}{the network's marginal log-likelihood}
  \item{N}{number of random networks to sample}
  \item{seed}{random seed}
  \item{type}{\code{mLL} or \code{FULLmLL} or \code{CONTmLL} or \code{CONTmLLDens} or \code{CONTmLLRatio}}
  \item{Pe}{prior of effect reporter positions in the phenotypic hierarchy (same dimension as D)}
  \item{para}{vector of length two: false positive rate and false negative rate for non-binary data. Used by \code{mLL}} 
  \item{hyperpara}{vector of length four: used by \code{FULLmLL()} for binary data}  
  \item{delta}{regularization parameter for automated E-gene subset selection (CONTmLLRatio only)}
  \item{selEGenes}{automated E-gene subset selection (includes tuning of delta for CONTmLLRatio)}  
}
\details{
Given data, N random network hypotheses from a null distribution are drawn as follows:  For each S-gene $S_k$ we randomly choose a number o of outgoing edges between 0 and 3. We then select o S-genes having at most 1 ingoing edge, connected  $S_k$ to them and transitively closed the graph. 
For all random network hypotheses it is counted, how often their likelihood is bigger than that of the given network. This yields an exact p-value.

Comparison of the likelihood of the given network to an empty and a fully connected gives us Bayes factors, which are also returned by the method.
}
\value{
\item{p.value}{p-value of the network according to the null hypothesis that it is random}
\item{BFzero}{marginal log-likelihood difference (weight of evidence) to the empty network}
\item{BFcomplete}{marginal log-likelihood difference (weight of evidence) to the fully connected network}
}

\author{Holger Froehlich}

\seealso{\code{\link{nem.bootstrap}}, \code{\link{nem}}}
\examples{
\dontrun{
   data("BoutrosRNAi2002")
   D <- BoutrosRNAiDiscrete[,9:16]
   p <- c(.13,.05)
   res = nem(D, para=p) # get best network
   nem.calcSignificance(D,as(res$graph,"matrix"),res$mLL, para=p) # assess its statistical significance
}
}

\keyword{models}
