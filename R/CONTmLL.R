# prob.inf - probability that Egene is influenced by specific Sgene
CONTmLL <- function (Phi, prob.inf, Pe) 
{
    if (!all(diag(Phi) == 1)) 
        stop("\nnem:mLL> Model main diagonal must be 1!")
        
    Phi <- Phi[unique(colnames(prob.inf)),unique(colnames(prob.inf))]
    nrS <- length(unique(colnames(prob.inf)))
    nrE <- nrow(prob.inf)
    nrA <- ncol(prob.inf)
        
    Sgenes <- unique(colnames(prob.inf))
    L <- matrix(nrow=nrE, ncol=nrS)
    
    # calculate likelihood
    for (i in seq(1,nrE)) {        
        for (j in seq(1,nrS)) {
            p <- 1
            for (k in seq(1,nrA)) {
                if (Phi[ceiling(k/2),j]==1)
                    p <- p * prob.inf[i,k]
                else
                    p <- p * (1-prob.inf[i,k])
            }
            L[i,j] <- p
        }
    }

    LP <- L * Pe
    s <- sum(log(rowSums(LP)))
    ep <- LP/rowSums(LP)
    colnames(ep) <- Sgenes
    rownames(ep) <- rownames(prob.inf)
    map <- Sgenes[apply(ep, 1, which.max)]
    names(map) <- rownames(prob.inf)
    
    return(list(mLL = s, pos = ep, mappos = map))
}
