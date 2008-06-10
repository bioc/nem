\name{score}
\alias{score}
\alias{print.score}
\alias{PhiDistr}

\title{Computes the marginal likelihood of phenotypic hierarchies}
\description{
   Function to compute the marginal likelihood of a set of phenotypic hierarchies.
}
\usage{
score(models, D, type="mLL", para=NULL, hyperpara=NULL, Pe=NULL, Pm=NULL, lambda=0, delta=1, verbose=TRUE, graphClass="graphNEL")

\method{print}{score}(x, ...)

PhiDistr(Phi, Pm, a=1, b=0.5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{models}{ a list of adjacency matrices with unit main diagonal }
  \item{D}{ data matrix. Columns correspond to the nodes in the
  silencing scheme. Rows are effect reporters. }
  \item{type}{\code{mLL} or \code{FULLmLL} or \code{CONTmLL} or \code{CONTmLLBayes} or \code{CONTmLLMAP}. \code{CONTmLLDens} and \code{CONTmLLRatio} are identical to \code{CONTmLLBayes} and \code{CONTmLLMAP} and are still supported for compatibility reasons. \code{mLL} and \code{FULLmLL} are used for binary data (see \code{BoutrosRNAiDiscrete}) and \code{CONTmLL} for a matrix of effect probabilities. \code{CONTmLLBayes} and \code{CONTmLLMAP} are used, if log-odds ratios, p-value densities or any other model specifies effect likelihoods. \code{CONTmLLBayes} refers to an inference scheme, were the linking positions of E-genes to S-Genes are integrated out, and \code{CONTmLLMAP} to an inference scheme, were a MAP estimate for the linking positions is calculated.}
  \item{para}{Vector with parameters \code{a} and \code{b} (for "mLL" with count data)}
  \item{hyperpara}{Vector with hyperparameters \code{a0}, \code{b0}, \code{a1}, \code{b1} for "FULLmLL"}
  \item{Pe}{prior position of effect reporters. Default: uniform over nodes in silencing scheme}
  \item{Pm}{prior on model graph (n x n matrix) with entries 0 <= priorPhi[i,j] <= 1 describing the probability of an edge between gene i and gene j.}
  \item{lambda}{regularization parameter to incorporate prior assumptions.}
  \item{delta}{regularization parameter for automated E-gene subset selection (CONTmLLRatio only)}
  \item{verbose}{output while running or not}
  \item{graphClass}{output inferred graph either as graphNEL or matrix}
  
  \item{x}{nem object}
  \item{...}{other arguments to pass}

  \item{Phi}{adjacency matrix}  
  \item{a}{parameter of the inverse gamma prior for v=1/lambda}
  \item{b}{parameter of the inverse gamma prior for v=1/lambda}
}
\details{
  
  Scoring models by marginal log-likelihood is implemented in function
  \code{score}. Input consists of models and data, the type of the score
  (\code{"mLL"}, \code{"FULLmLL"}, \code{"CONTmLL"} or \code{"CONTmLLBayes"} or \code{"CONTmLLMAP"}), the corresponding paramters (\code{para}) or hyperparameters (\code{hyperpara}), a prior for phenotype
  positions (\code{Pe}) and model structures \code{Pm} with regularization parameter \code{lambda}. If a structure prior \code{Pm} is provided, but no regularization parameter \code{lambda}, Bayesian model averaging with an inverse gamma prior on 1/lambda is performed. 
  With type "CONTmLLMAP" usually an automated selection of most relevant E-genes is performed by introducing a "null" S-gene. The corresponding prior probability of leaving out an E-gene is set to delta/no. S-genes. 
  
  \code{score} is usually called within function \code{nem}.

}
\value{
  nem object
}
\references{
\item{[1]}{Markowetz F, Bloch J, Spang R, Non-transcriptional pathway features reconstructed from secondary effects of RNA interference, Bioinformatics, 2005.}
\item{[2]}{Markowetz F, Probabilistic Models for Gene Silencing Data, PhD thesis, Free University Berlin, 2006.}
\item{[3]}{Froehlich H, Fellmann M, Sueltmann H, Poustka A, Beissbarth T: Estimating Large Scale Signaling Networks through Nested Effects Models from Intervention Effects in Microarray Data. Bioinformatics, 1, 2008.}
\item{[4]}{Froehlich H, Fellmann M, Sueltmann H, Poustka A, Beissbarth T: Large Scale Statistical Inference of Signaling Pathways from RNAi and Microarray Data, BMC Bioinformatics, 8:386, 2007.}
}
\author{Holger Froehlich <URL: http://www.dkfz.de/mga2/people/froehlich>, Florian Markowetz <URL: http://genomics.princeton.edu/~florian>}

\seealso{\code{\link{nem}}, \code{\link{mLL}}, \code{\link{FULLmLL}}, \code{\link{enumerate.models}} }
\examples{
   # Drosophila RNAi and Microarray Data from Boutros et al, 2002
   data("BoutrosRNAi2002")
   D <- BoutrosRNAiDiscrete[,9:16]

   # enumerate all possible models for 4 genes
   models <- enumerate.models(unique(colnames(D)))

   # score models with marginal likelihood
   result <- score(models,D,type="mLL",para=c(.13,.05))
   
   # plot graph of the best model
   plot(result,what="graph")

   # plot scores
   plot(result,what="mLL") 
    
   # plot posterior of E-gene positions according to best model
   plot(result,what="pos")
   
   # MAP estimate of effect positions for the best model
   result$mappos

}

\keyword{models}
