\name{plot.nem}
\alias{plot.nem}
\alias{plot.triples}
\alias{plot.pairwise}
\alias{plot.score}
\alias{plot.ModuleNetwork}

\title{plot nested effect model}
\description{
  plot graph of nested effect model, the marginal likelihood distribution or the posterior position of the effected genes
}
\usage{
	plot.nem(x, what="graph", remove.singletons=FALSE, PDF=FALSE, filename="nemplot.pdf", thresh=0, transitiveReduction=FALSE, ...)
}

\arguments{
 	\item{x}{nem object to plot}
	\item{what}{(i), "graph", (ii) "mLL" = likelihood distribution, (iii) "pos" = posterior position of effected genes}
	\item{remove.singletons}{remove unconnected nodes from the graph plot}
	\item{PDF}{output as PDF-file}
	\item{filename}{filename of PDF-file}
	\item{thresh}{if x has a real valued adjacency matrix (weight matrix),  don't plot edges with weight <= thresh}
	\item{transitiveReduction}{plot a transitively reduced graph} 
	\item{...}{other arguments to pass}
}
\value{
  none
}


\author{Florian Markowetz <URL: http://genomics.princeton.edu/~florian>}

\seealso{\code{\link{nem}}, \code{\link{score}}}

\keyword{models}
