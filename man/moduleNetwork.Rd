\name{moduleNetwork}
\alias{moduleNetwork}
\alias{print.ModuleNetwork}

\title{Infers a phenotypic hierarchy using the module network}
\description{
  Function \code{moduleNetwork} estimates the hierarchy using a divide and conquer approach. In each step only a subset of nodes (called module)
  is involved and no exhaustive enumeration of model space is needed as in function \code{score}.
}
\usage{
moduleNetwork(D,type="mLL",Pe=NULL,Pm=NULL,lambda=0,para=NULL,hyperpara=NULL,selEGenes=FALSE,verbose=TRUE)

#S3 methods for class 'moduleNetwork'
print.ModuleNetwork(x,...)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{D}{data matrix. Columns correspond to the nodes in the silencing scheme. Rows are phenotypes.}   
\item{type}{(1.) marginal likelihood "mLL" (only for cout matrix D), or (2.) full marginal likelihood "FULLmLL" integrated over a and b and depending on hyperparameters a0, a1, b0, b1 (only for count matrix D), or (3.) "CONTmLL" marginal likelihood for probability matrices, or (4.) "CONTmLLDens" marginal likelihood for probability density matrices, or (5.) "CONTmLLRatio" for log-odds ratio matrices}  
  \item{Pe}{prior position of effect reporters. Default: uniform over nodes in hierarchy}  
  \item{Pm}{prior on model graph (n x n matrix) with entries 0 <= priorPhi[i,j] <= 1 describing the probability of an edge between gene i and gene j.}
  \item{lambda}{regularization parameter to incorporate prior assumptions.}  
  \item{para}{vector with parameters a and b for "mLL", if count matrices are used}
  \item{hyperpara}{vector with hyperparameters a0, b0, a1, b1 for "FULLmLL"}  
  \item{selEGenes}{optimize selection of E-genes for each model}
  \item{verbose}{do you want to see progress statements printed or not? Default: TRUE}  

  \item{x}{nem object}
  \item{...}{other arguments to pass}
}

\details{
  \code{moduleNetwork} is an alternative to exhaustive search 
  by the function \code{score} and more accurate than \code{pairwise.posterior}. 
  It uses clustering to sucessively split the network into smaller modules, which can then be estimated completely. Connections between modules are estimated pairwise between nodes of two different modules.
}

\value{
  \item{graph}{the inferred directed graph (graphNEL object)}    
  \item{pos}{posterior over effect positions}
  \item{mappos}{MAP estimate of effect positions} 
  \item{type}{as used in function call}
  \item{para}{as used in function call}
  \item{hyperpara}{as used in function call}
  \item{lambda}{as in function call}
}

\references{
  Froehlich H, Fellmann M, Sueltmann H, Poustka A, Beissbarth T:
  Large Scale Statistical Inference of Signaling Pathways from RNAi and
  Microarray Data.
  BMC Bioinformatics, 2007.

  Froehlich H, Fellmann M, Sueltmann H, Poustka A, Beissbarth T:
  Estimating Large Scale Signaling Networks through Nested Effects
  Models from Intervention Effects in Microarray Data.
  Proc. German Conf. Bioinformatics (GCB), pp. 45 - 54, 2007.  
}
  
\author{Holger Froehlich}


\seealso{\code{\link{score}}, \code{\link{nem}}}
\examples{
   data("BoutrosRNAi2002") 
   res <- moduleNetwork(BoutrosRNAiDiscrete[,9:16],para=c(.13,.05))
   
   # plot graph
   plot(res,what="graph")
   
   # plot posterior over effect positions
   plot(res,what="pos")
   
   # estimate of effect positions
   res$mappos
   
}
\keyword{graphs}
