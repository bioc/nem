\name{score}
\alias{score}
\alias{plot.score}
\alias{print.score}

\title{Computes the marginal likelihood of phenotypic hierarchies}
\description{
   Function to compute the marginal likelihood of a set of phenotypic hierarchies.
}
\usage{
score(models, D, type = "mLL", para = NULL, hyperpara = NULL, Pe = NULL, verbose = TRUE)

# S3 metehods for class 'score'
plot.score(x, what="graph",remove.singletons=FALSE, PDF=FALSE, filename="nemplot.pdf", ...)
print.score(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{models}{ a list of adjacency matrices with unit main diagonal }
  \item{D}{ data matrix. Columns correspond to the nodes in the
  silencing scheme. Rows are effect reporters. }
  \item{type}{Marginal likelihood "mLL" depending on parameters \code{a} and \code{b}, or
  full marginal likelihood "FULLmLL" integrated over \code{a} and \code{b} and depending on
  hyperparameters \code{a0}, \code{a1}, \code{b0}, \code{b1}. }
  \item{para}{Vector with parameters \code{a} and \code{b} }
  \item{hyperpara}{Vector with hyperparameters \code{a0}, \code{b0}, \code{a1}, \code{b1}}
  \item{Pe}{Prior position of effect reporters. Default: uniform over nodes in silencing scheme}
  \item{verbose}{output while running or not}
  \item{x}{an object of class 'score'}
  \item{what}{type of plot: 'graph', 'mLL', or 'pos'. Default: 'graph'}
  \item{remove.singletons}{remove single nodes which are not connected to any other node when plotting? Default: FALSE}
  \item{PDF}{output as pdf file? Default: FALSE}
  \item{filename}{name of the pdf if any. Default: "nemplot.pdf"}
  \item{...}{additional arguments for plotting}
}
\details{
  
  Scoring models by marginal log-likelihood is implemented in function
  \code{score}. Input consists of models and data, the type of the score
  (\code{"mLL"} or \code{"FULLmLL"}), the corresponding paramters
  (\code{para}) or hyperparameters (\code{hyperpara}) and a prior for phenotype
  positions (\code{Pe}).
  
  \code{score} is usually called from within function \code{nem}.

}
\value{
  \item{graph}{the model with highest marginal likelihood (graphNEL object)}
  \item{mLL}{vector of marginal likelihoods for all models}
  \item{pos}{a list of estimated positions of effect reporters for each model}
  \item{mappos}{a list of maximum aposteriori estimates of effect positions for each model}
  \item{type}{as used in function call}
  \item{para}{as used in function call}
  \item{hyperpara}{as used in function call}
}
\references{
\item{[1]}{Markowetz F, Bloch J, Spang R, Non-transcriptional pathway features reconstructed from secondary effects of RNA interference, Bioinformatics, 2005}
\item{[2]}{Markowetz F, Probabilistic Models for Gene Silencing Data, PhD thesis, Free University Berlin, 2006}
}
\author{Florian Markowetz <URL: http://genomics.princeton.edu/~florian>}

\seealso{\code{\link{nem}}, \code{\link{mLL}}, \code{\link{FULLmLL}}, \code{\link{enumerate.models}} }
\examples{
   # Drosophila RNAi and Microarray Data from Boutros et al, 2002
   data("BoutrosRNAi2002")
   D <- BoutrosRNAiDiscrete[,9:16]

   # enumerate all possible models for 4 genes
   models <- enumerate.models(unique(colnames(D)))

   # score models with marginal likelihood
   result <- score(models,D,type="mLL",para=c(.13,.05))
   
   # plot graph
   plot(result,what="graph")

   # plot scores
   plot(result,what="mLL") 
    
   # plot posterior of E-gene positions
   plot(result,what="pos")
   
   # MAP estimate of effect positions
   result$mappos[[which.max(result$mLL)]]

}

\keyword{models}
