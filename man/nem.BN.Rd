\name{nem.BN}
\alias{nem.BN}
\alias{print.nem.BN}

\title{Bayesian Network Nested Effects Models}
\description{Uses a Bayesian network interpretation of Nested Effects Models to estimate the signals graph.}
\usage{
nem.BN(D, inference="greedy", mode="binary_ML", lambda=0, verbose=TRUE)
}

\arguments{
  \item{D}{data matrix with experiments in the columns (binary or continious)}
  \item{inference}{\code{exhaustive} to use exhaustive enumeration; or \code{greedy} for optimizing the linking of effects to signals and the signals graph in an alternating fashion}
  \item{mode}{\code{binary_ML}: effects come from a binomial distribution - ML learning of parameters; \code{binary_Bayesian}: effects come from a binomial distribution - Bayesian learning of parameters with beta distribution prior; \code{continous_ML}: effects come from a normal distribution - ML learning of parameters; \code{continous_Bayesian}: effects come from a normal distribution - Bayesian learning of parameters with gamma distribution prior.}  
  \item{lambda}{regularization parameter to incorporate prior assumptions. Not used so far.}  
  \item{verbose}{do you want to see progression statements" Default: TRUE}
}
\details{
  \code{plot.nem} plots the inferred phenotypic hierarchy as a directed graph.    
}
\value{
An object of class 'nem.BN'
\item{graph}{the inferred phenotypic hierarchy}
\item{mLL}{log (posterior) marginal likelihood}
\item{mappos}{estimated position of effects in the phenotypic hierarchy}
\item{selected}{selected E-gene subset}
\item{type}{= mode in function call}
\item{lambda}{see above}
}

\author{Cordula Zeller, Holger Froehlich <URL: http:/www.dkfz.de/mga2/people/froehlich>}


\seealso{\code{plot.nem}}
\keyword{graphs}
\keyword{models}
