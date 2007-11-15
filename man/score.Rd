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

PhiDistr(Phi, Pm, a=1, b=ncol(Phi)^2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{models}{ a list of adjacency matrices with unit main diagonal }
  \item{D}{ data matrix. Columns correspond to the nodes in the
  silencing scheme. Rows are effect reporters. }
  \item{type}{(1.) marginal likelihood "mLL" (only for cout matrix D), or (2.) full marginal likelihood "FULLmLL" integrated over a and b and depending on hyperparameters a0, a1, b0, b1 (only for count matrix D), or (3.) "CONTmLL" marginal likelihood for probability matrices, or (4.) "CONTmLLDens" marginal likelihood for probability log-density matrices, or (5.) "CONTmLLRatio" for log-odds ratio matrices}  
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
  (\code{"mLL"}, \code{"FULLmLL"}, \code{"CONTmLL"} or \code{"CONTmLLDens"} or \code{"CONTmLLRatio"}), the corresponding paramters (\code{para}) or hyperparameters (\code{hyperpara}), a prior for phenotype
  positions (\code{Pe}) and model structures \code{Pm} with regularization parameter \code{lambda}. If a structure prior \code{Pm} is provided, but no regularization parameter \code{lambda}, Bayesian model averaging with an inverse gamma prior on 1/lambda is performed. 
  With type "CONTmLLRatio" usually an automated selection of most relevant E-genes is performed by introducing a "null" S-gene. The corresponding prior probability of leaving out an E-gene is set to delta/no. S-genes. 
  
  \code{score} is usually called within function \code{nem}.

}
\value{
  \item{graph}{the inferred directed graph (graphNEL object)}   
  \item{mLL}{marginal likelihood of final model} 
  \item{pos}{posterior over effect positions}
  \item{mappos}{MAP estimate of effect positions} 
  \item{type}{as used in function call}
  \item{para}{as used in function call}
  \item{hyperpara}{as used in function call}
  \item{lambda}{as in function call}
  \item{selected}{selected E-gene subset}
}
\references{
\item{[1]}{Markowetz F, Bloch J, Spang R, Non-transcriptional pathway features reconstructed from secondary effects of RNA interference, Bioinformatics, 2005.}
\item{[2]}{Markowetz F, Probabilistic Models for Gene Silencing Data, PhD thesis, Free University Berlin, 2006.}
\item{[3]}{Froehlich H, Fellmann M, Sueltmann H, Poustka A, Beissbarth T: Estimating Large Scale Signaling Networks through Nested Effects Models from Intervention Effects in Microarray Data. Proc. German Conf. Bioinformatics (GCB), pp. 45 - 54, 2007.}
\item{[4]}{Froehlich H, Fellmann M, Sueltmann H, Poustka A, Beissbarth T: Large Scale Statistical Inference of Signaling Pathways from RNAi and Microarray Data, BMC Bioinformatics, 8:386, 2007.}
}
\author{Holger Froehlich, Florian Markowetz <URL: http://genomics.princeton.edu/~florian>}

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
