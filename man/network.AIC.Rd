\name{network.AIC}
\alias{network.AIC}

\title{AIC criterion for network graph}
\description{
  calclate AIC for a given network graph (should be transitively closed)
}
\usage{
	network.AIC(network,k=2,verbose=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{network}{a nem object (e.g. 'pairwise')}
  \item{k}{penalty per parameter in the AIC calculation. Default = 2 for classical AIC}
  \item{verbose}{print out the result}  
}

\details{
	For k = log(n) the BIC (Schwarz criterion) is computed. Usually this function is not called directly but from \code{nemModelSelection}
}

\value{
	AIC value
 
}

\author{Holger Froehlich}


\seealso{\code{\link{nemModelSelection}}}
\examples{
   data("BoutrosRNAi2002") 
   res1 <- moduleNetwork(BoutrosRNAiDiscrete[,9:16],para=c(.13,.05))
   network.AIC(res1)
   res2 <- moduleNetwork(BoutrosRNAiDiscrete[,9:16],para=c(.13,.05),Pm=matrix(0,ncol=4,nrow=4),lambda=10)
   network.AIC(res2)
}
\keyword{graphs}
