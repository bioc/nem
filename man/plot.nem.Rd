\name{plot.nem}
\alias{plot.nem}
\alias{plot.triples}
\alias{plot.pairwise}
\alias{plot.score}
\alias{plot.ModuleNetwork}
\alias{plot.nem.greedy}
\alias{plot.nem.greedyMAP}
\alias{plot.nem.bootstrap}
\alias{plot.nem.jackknife}
\alias{plot.nem.consensus}

\title{plot nested effect model}
\description{
  plot graph of nested effect model, the marginal likelihood distribution or the posterior position of the effected genes
}
\usage{
	\method{plot}{nem}(x, what="graph", remove.singletons=FALSE, PDF=FALSE, filename="nemplot.pdf", thresh=0, transitiveReduction=FALSE, plot.probs=FALSE, SCC=TRUE, ...)
}

\arguments{
 	\item{x}{nem object to plot}
	\item{what}{(i), "graph", (ii) "mLL" = likelihood distribution, (iii) "pos" = posterior position of effected genes}
	\item{remove.singletons}{remove unconnected nodes from the graph plot}
	\item{PDF}{output as PDF-file}
	\item{filename}{filename of PDF-file}
	\item{thresh}{if x has a real valued adjacency matrix (weight matrix),  don't plot edges with |weight| <= thresh}
	\item{transitiveReduction}{plot a transitively reduced graph} 
 	\item{plot.probs}{plot edge weights/probabilities. If regulation directions have been inferred (see \code{infer.edge.type}), upregulated edges are drawn blue and downregulated edges red. Edges, were no clear direction could be inferred, are drawn in black.} 
	\item{SCC}{plot the strongly connected components graph}
	\item{...}{other arguments to pass}
}
\value{
  none
}


\author{Florian Markowetz <URL: http://genomics.princeton.edu/~florian>, Holger Froehlich <URL: http://www.dkfz.de/mga2/people/froehlich>}

\seealso{\code{\link{nem}}, \code{\link{plot.effects}}, \code{\link{infer.edge.type}}}

\keyword{models}
