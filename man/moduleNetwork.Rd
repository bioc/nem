\name{moduleNetwork}
\alias{moduleNetwork}
\alias{print.ModuleNetwork}

\title{Infers a phenotypic hierarchy using the module network}
\description{
  Function \code{moduleNetwork} estimates the hierarchy using a divide and conquer approach. In each step only a subset of nodes (called module)
  is involved and no exhaustive enumeration of model space is needed as in function \code{score}.
}
\usage{
moduleNetwork(D,type="mLL",Pe=NULL,Pm=NULL,lambda=0,delta=1,para=NULL,hyperpara=NULL,verbose=TRUE)

\method{print}{ModuleNetwork}(x,...)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{D}{data matrix. Columns correspond to the nodes in the silencing scheme. Rows are phenotypes.}   
  \item{type}{see \code{nem}}
  \item{Pe}{prior position of effect reporters. Default: uniform over nodes in hierarchy}  
  \item{Pm}{prior on model graph (n x n matrix) with entries 0 <= priorPhi[i,j] <= 1 describing the probability of an edge between gene i and gene j.}
  \item{lambda}{regularization parameter to incorporate prior assumptions.}  
  \item{delta}{regularization parameter for automated E-gene subset selection (CONTmLLRatio only)}
  \item{para}{vector with parameters a and b for "mLL", if count matrices are used}
  \item{hyperpara}{vector with hyperparameters a0, b0, a1, b1 for "FULLmLL"}    
  \item{verbose}{do you want to see progress statements printed or not? Default: TRUE}  

  \item{x}{nem object}
  \item{...}{other arguments to pass}
}

\details{
  \code{moduleNetwork} is an alternative to exhaustive search 
  by the function \code{score} and more accurate than \code{pairwise.posterior} and \code{triples.posterior}. 
  It uses clustering to sucessively split the network into smaller modules, which can then be estimated completely. Connections between modules are estimated by performing a constraint greedy hillclimbing.
}

\value{
 nem object
}

\references{
  \item{[1]}{Froehlich H, Fellmann M, Sueltmann H, Poustka A, Beissbarth T:
  Large Scale Statistical Inference of Signaling Pathways from RNAi and
  Microarray Data.
  BMC Bioinformatics, 2007.}
  \item{[2]}{Froehlich H, Fellmann M, Sueltmann H, Poustka A, Beissbarth T:
  Estimating Large Scale Signaling Networks through Nested Effects
  Models from Intervention Effects in Microarray Data.
  Bioinformatics, 1, 2008. }
}
  
\author{Holger Froehlich}


\seealso{\code{\link{score}}, \code{\link{nem}}}
\examples{
   data("BoutrosRNAi2002") 
   res <- nem(BoutrosRNAiDiscrete[,9:16],para=c(.13,.05),inference="ModuleNetwork")
   
   # plot graph
   plot(res,what="graph")
   
   # plot posterior over effect positions
   plot(res,what="pos")
   
   # estimate of effect positions
   res$mappos
   
}
\keyword{graphs}
