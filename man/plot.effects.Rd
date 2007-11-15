\name{plot.effects}
\alias{plot.effects}

\title{Plots data according to a phenotypic hierarchy}
\description{
  \code{plot.effects} visualizes the subset structure in the data by reordering rows and columns according to the topological order given by a phenotypic hierarchy.
}
\usage{
\method{plot}{effects}(x,nem,border=TRUE,...)
}

\arguments{
  \item{x}{data matrix}
  \item{nem}{phenotypic hierarchy (object of class 'score' or 'pairwise')}
  \item{border}{draw red lines to indicate gene-specific effect reporters. Default: TRUE}
  \item{...}{additional parameters for the plotting function}
}
\details{
  The experiments in the columns are reordered according to the topological order given by a phenotypic hierarchy.
  The effect reporters in the rows are grouped together by their position in the hierarchy.
  The groups are then arranged by topological order.
  Within each group the rows are hierarchically clustered.
  }

\author{Florian Markowetz <URL: http://genomics.princeton.edu/~florian>}

\examples{
   data("BoutrosRNAi2002")
   D <- BoutrosRNAiDiscrete[,9:16]
   res <- nem(D,para=c(.13,.05))
   plot.effects(D,res)
}
\keyword{graphs}
