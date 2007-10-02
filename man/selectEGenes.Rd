\name{getRelevantEGenes}
\alias{getRelevantEGenes}
\alias{selectEGenes}
\alias{filterEGenes}

\title{Automatic selection of most relevant S-genes}
\description{
   1. Selects those E-genes, which have the highest likelihood under the given network hypothesis.
   2. Cluster E-genes and select one E-gene from each cluster to reduce the amount of data on the array.
}
\usage{
filterEGenes(Porig, D, ntop=100)

getRelevantEGenes(Phi, D, nEgenes=min(5*ncol(Phi), nrow(D1)), type="mLL", para=NULL, hyperpara=NULL, Pe=NULL, Pm=NULL, lambda=0)

selectEGenes(Phi,D1,D0=NULL,para=NULL,hyperpara=NULL,Pe=NULL,Pm=NULL,lambda=0,type="mLL", nEgenes=min(5*ncol(Phi), nrow(D1)))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Porig}{matrix of raw p-values, typically from the complete array}
  \item{D}{ data matrix. Columns correspond to the nodes in the silencing scheme. Rows are effect reporters. }
  \item{ntop}{number of top genes to consider from each knock-down experiment}

  \item{Phi}{ adjacency matrix with unit main diagonal }    
  \item{nEgenes}{ no. of E-genes to select}
  \item{type}{(1.) marginal likelihood "mLL" (only for cout matrix D), or (2.) full marginal likelihood "FULLmLL" integrated over a and b and depending on hyperparameters a0, a1, b0, b1 (only for count matrix D), or (3.) "CONTmLL" marginal likelihood for probability matrices, or (4.) "CONTmLLDens" marginal likelihood for probability density matrices, or (5.) "CONTmLLRatio" for log-odds ratio matrices}  
  \item{para}{Vector with parameters \code{a} and \code{b} (for "mLL" with count data)}
  \item{hyperpara}{Vector with hyperparameters \code{a0}, \code{b0}, \code{a1}, \code{b1} for "FULLmLL"}
  \item{Pe}{prior position of effect reporters. Default: uniform over nodes in silencing scheme}
  \item{Pm}{prior on model graph (n x n matrix) with entries 0 <= priorPhi[i,j] <= 1 describing the probability of an edge between gene i and gene j.}
  \item{lambda}{regularization parameter to incorporate prior assumptions.}
  
   \item{D1}{(i) count matrix for discrete data: phenotypes x genes. How often did we see an effect after interventions? (ii) matrix describing the probabilities of an effect (iii) probability density matrix discribing the strength of an effect}
  \item{D0}{count matrix: phenotypes x genes. How often did we NOT see an effect after intervention? Not used for continious data}
}
\details{
	uses \code{mLL} or \code{FULLmLL} to score each E-gene.
}
\value{
  \item{I}{index of selected E-genes}

  \item{mLL}{marginal likelihood of a phenotypic hierarchy}
  \item{pos}{posterior distribution of effect positions in the hierarchy}
  \item{mappos}{Maximum aposteriori estimate of effect positions}
}

\author{Holger Froehlich}

\seealso{\code{\link{nem}}, \code{score}, \code{\link{mLL}}, \code{\link{FULLmLL}}, \code{\link{enumerate.models}} }
\examples{
   # Drosophila RNAi and Microarray Data from Boutros et al, 2002
   data("BoutrosRNAi2002")
   D <- BoutrosRNAiDiscrete[,9:16]

   # enumerate all possible models for 4 genes
   models <- enumerate.models(unique(colnames(D)))  
   
   getRelevantEGenes(models[[64]], D, para=c(.13,.05))

}

\keyword{models}
