createBN = function(data, datatotal=NULL, poss.vertices = NULL,coregraph = NULL, nullnode  = TRUE, marginal = NULL){
	if (!is.character(colnames(data))) {
        colnames(data) = as.character(1:ncol(data))
    }
	if (!is.character(rownames(data))) {
        rownames(data) = as.character(1:nrow(data))
    }
	if(is.null(datatotal)){
		datatotal = matrix(1, nrow = nrow(data), ncol = ncol(data))
	}
	if(is.null(poss.vertices)){
		poss.vertices = 1:length(unique(colnames(data)))
	}
	dimnames(datatotal) = dimnames(data)
	cores = unique(colnames(data))
	if(nullnode){
		cores = c(cores, "null")
		}
	reporters = rownames(data)
	nrrep = length(reporters)
	nrcores = length(cores)
	marginals = rownames(data)
	if(is.null(coregraph)){
		coregraph = matrix(0,ncol = nrcores, nrow = nrcores)}
	if(is.null(marginal)){
		marginal = matrix(0,ncol = nrrep, nrow = nrcores)
	}
	dimnames(coregraph) = list(cores, cores)
	dimnames(marginal) = list(cores, marginals)
	diag(coregraph) = 1
	if(nullnode){coregraph["null",]=1
		coregraph[,"null"]=1
		}
	wh = which(coregraph == 0)
	if(nullnode){coregraph["null",] = 0
		coregraph[,"null"]=0
		}
	parameters = matrix(1,nrow=4,ncol=nrrep)
	ve = 1
	k=1
	for(i in wh){
		row = i %% nrow(coregraph)
		if (row==0){row = nrow(coregraph)}
		col = ceiling(i/nrow(coregraph))
		common = (data[,row]%*% data[,col])
		if (common > 0){ve[k]=i
			k = k + 1}}
	BN = list(data = data, datatotal = datatotal, coregraph = coregraph, marginal = marginal, parameters = parameters, 
			poss.vertices = poss.vertices, variable.edges = ve, score = -Inf)
	return(BN)
}

fit.BN = function(BN, method = "greedy", verbose = FALSE, trace = FALSE, mode = "binary_ML",...){
	switch(method, exhaustive = {
            return(exhaustive_BN(BN, verbose = verbose, trace = trace, 
            mode = mode,...))
        }, greedy = {
            return(ingreed_BN(BN, verbose = verbose, trace = trace,
               mode = mode,...))
        }, stop("This method does not exist"))
}

which.is.max = function (x){
    y <- seq_along(x)[x == max(x)]
    if (length(y) > 1)
        sample(y, 1)
    else y
}

# Vor Aufruf der score Funktionen muss gesichert sein, dass der coregraph transitiv abgeschlossen ist und dass die parameters Parameter aktualisiert sind. 
score_discrete_ML = function(BN, loc = NULL){
	parameters = BN$parameters
	if(!is.null(loc)){
		wholegraph = ((BN$coregraph %*% BN$marginal[,loc]) > 0) * 1
		names(wholegraph) = rownames(BN$marginal)
		wholegraph = wholegraph[colnames(BN$data)]
		score1 = dbinom(sum(BN$data[loc, ] * wholegraph), sum(BN$datatotal[loc, ] * wholegraph), parameters[4, loc], log = TRUE)
		score0 = dbinom(sum(BN$data[loc, ] * (1 - wholegraph)), sum(BN$datatotal[loc,] * (1 - wholegraph)), parameters[2,loc], log = TRUE)
		return(score1 + score0)
	}
	else{
		wholegraph = ((BN$coregraph %*% BN$marginal) > 0) * 1
		rownames(wholegraph) = rownames(BN$coregraph)
		wholegraph = wholegraph[colnames(BN$data),]
		score = matrix(1, nrow = nrow(BN$data), ncol = 2)
		for(reporter in 1:nrow(BN$data)){
			score[reporter, 1] = dbinom(BN$data[reporter,] %*% wholegraph[,reporter],BN$datatotal[reporter,] %*% wholegraph[, reporter], parameters[4,reporter], log = TRUE)
			score[reporter, 2] = dbinom(BN$data[reporter,] %*% (1 - wholegraph[,reporter]),BN$datatotal[reporter,] %*% (1 - wholegraph[,reporter]), parameters[2,reporter], log = TRUE)
		}
		sc = sum(score)
		return(sum(score))}
}

parameters_discrete_ML = function(BN, loc = NULL){
	if(!is.null(loc)){
		parameters = numeric(4)
		wholegraph = ((BN$coregraph %*% BN$marginal[,loc]) > 0) * 1
		names(wholegraph) = rownames(BN$marginal)
		wholegraph = wholegraph[colnames(BN$data)]
		predicted.pos = BN$datatotal[loc,] %*% wholegraph		
		if (predicted.pos == 0){predicted.pos = 1}
		parameters[4] = max(0.5, (BN$data[loc,] %*% wholegraph) / predicted.pos)
		parameters[3] = 1 - parameters[4]
		predicted.neg = BN$datatotal[loc, ] %*% (1 - wholegraph)
		if (predicted.neg == 0){predicted.neg = 1}
		parameters[1] = max(0.5, ((BN$datatotal[loc,] - BN$data[loc,]) %*% (1 - wholegraph)) / predicted.neg)
		parameters[2] = 1 - parameters[1]
		return(parameters)
	}
	else{parameters = BN$parameters
		wholegraph = ((BN$coregraph %*% BN$marginal) > 0) * 1
		rownames(wholegraph) = rownames(BN$coregraph)
		wholegraph = wholegraph[colnames(BN$data),]
		score = matrix(1, nrow = nrow(BN$data), ncol = 2)
		for(reporter in 1:nrow(BN$data)){
			predicted.pos = BN$datatotal[reporter,] %*% wholegraph[,reporter]		
			if (predicted.pos == 0){predicted.pos = 1}
			parameters[4,reporter] = max(0.5,(BN$data[reporter,]%*%wholegraph[,reporter])/predicted.pos)
			parameters[3,reporter] = 1 - parameters[4,reporter]
			predicted.neg = BN$datatotal[reporter, ] %*% (1 - wholegraph[,reporter])
			if (predicted.neg == 0){predicted.neg = 1}	
			parameters[1, reporter] = max(0.5,((BN$datatotal[reporter,] - BN$data[reporter,]) %*% (1 - wholegraph[,reporter])) / predicted.neg)
			parameters[2, reporter] = 1 - parameters[1,reporter]
		}
		colnames(parameters) = rownames(BN$data)
		rownames(parameters) = c("p00","p01","p10","p11")
		return(parameters)}
}

parameters_discrete_Bayesian=function(BN, alpha1 = 5, beta1 = 2, alpha0 = 2, beta0 = 5){
		parameters = BN$parameters
		wholegraph = ((BN$coregraph %*% BN$marginal) > 0) * 1
		rownames(wholegraph) = rownames(BN$coregraph)
		wholegraph = wholegraph[colnames(BN$data),]
		for(reporter in 1:nrow(BN$data)){
			predicted.pos = BN$datatotal[reporter,] %*% wholegraph[,reporter]		
			parameters[4, reporter] = alpha1 + BN$data[reporter,] %*% wholegraph[,reporter]
			parameters[3, reporter] = beta1 + predicted.pos - parameters[4,reporter]
			predicted.neg = BN$datatotal[reporter,] %*% (1 - wholegraph[,reporter])		
			parameters[1, reporter] = beta0 +  (BN$datatotal[reporter,] - BN$data[reporter,]) %*% (1 - wholegraph[,reporter])
			parameters[2, reporter] = alpha0 + predicted.neg - parameters[1,reporter]
		}
		colnames(parameters) = rownames(BN$data)
		rownames(parameters) = c("beta0", "alpha0", "beta1", "alpha1")
		return(parameters)
	}

score_discrete_Bayesian = function(BN, loc = NULL, alpha1 = 5, beta1 = 2, alpha0 = 2, beta0 = 5, likelihood = FALSE){
	if(!is.null(loc)){
		parameters = numeric(4)
		wholegraph = ((BN$coregraph %*% BN$marginal[,loc]) > 0) * 1
		names(wholegraph) = rownames(BN$marginal)
		wholegraph = wholegraph[colnames(BN$data)]
		predicted.pos = BN$datatotal[loc, ] %*% wholegraph
		parameters[4] = BN$data[loc,] %*% wholegraph
		parameters[3] = predicted.pos - parameters[4]
		predicted.neg = BN$datatotal[loc, ] %*% (1 - wholegraph)
		parameters[1] = (BN$datatotal[loc,] - BN$data[loc,]) %*% (1 - wholegraph)
		parameters[2] = predicted.neg - parameters[1]	
		score = lfactorial(parameters[4] + alpha1 - 1) + lfactorial(parameters[3] + beta1 - 1) + lfactorial(parameters[2] + alpha0 - 1) + lfactorial(parameters[1] + beta0 - 1) - (
						lfactorial(predicted.pos + alpha1 + beta1 - 1) + lfactorial(predicted.neg + alpha0 + beta0 - 1))
		return(score)
	}
	else{
		parameters = BN$parameters
		wholegraph = ((BN$coregraph %*% BN$marginal) > 0) * 1
		rownames(wholegraph) = rownames(BN$coregraph)
		wholegraph = wholegraph[colnames(BN$data),]
		score = matrix(1, ncol = 1, nrow = nrow(BN$data))
		for(reporter in 1:nrow(BN$data)){
			predicted.pos = BN$datatotal[reporter,] %*% wholegraph[,reporter]		
			parameters[4,reporter] = BN$data[reporter,] %*% wholegraph[,reporter]
			parameters[3,reporter] = predicted.pos - parameters[4,reporter]
			predicted.neg = BN$datatotal[reporter,] %*% (1 - wholegraph[,reporter])		
			parameters[1,reporter] = (BN$datatotal[reporter,] - BN$data[reporter,]) %*% (1 - wholegraph[,reporter])
			parameters[2,reporter] = predicted.neg - parameters[1,reporter]
			score[reporter] = lfactorial(parameters[4, reporter] + alpha1 - 1) + lfactorial(parameters[3, reporter] + beta1 - 1) + lfactorial(parameters[2, reporter] + alpha0 - 1) + lfactorial(parameters[1,reporter] + beta0 - 1) - (
							lfactorial(predicted.pos + alpha1 + beta1 - 1) + lfactorial(predicted.neg + alpha0 + beta0 - 1))
		
		}
		if(likelihood){
			score = score + lfactorial(alpha0 + beta0 - 1) + lfactorial(alpha1 + beta1 - 1) - lfactorial(alpha0 - 1)- lfactorial(alpha1 - 1) - lfactorial(beta0 - 1)- lfactorial(beta1 - 1)
			}
		return(sum(score))
	}
}

parameters_continuous_ML = function(BN, loc = NULL){
	parameters = BN$parameters
	if(!is.null(loc)){
		parameters = numeric(4)
		wholegraph = ((BN$coregraph %*% BN$marginal[,loc]) > 0) * 1
		names(wholegraph) = rownames(BN$marginal)
		wholegraph = wholegraph[colnames(BN$data)]
		predicted.pos = sum(wholegraph)
		if (predicted.pos==0){predicted.pos=1}
		parameters[3] = (BN$data[loc,] %*% wholegraph)/predicted.pos
		parameters[4] = sum(((BN$data[loc,] - parameters[3])^2) %*% wholegraph) / predicted.pos
		predicted.neg = sum(1 - wholegraph)
		if(predicted.neg == 0){predicted.neg = 1}
		parameters[1] = (BN$data[loc,]%*%(1 - wholegraph)) / predicted.neg
		parameters[2] = sum(((BN$data[loc,] - parameters[1])^2) %*% (1 - wholegraph)) / predicted.neg
		return(parameters)
	}
	else{
		parameters = BN$parameters
		wholegraph = ((BN$coregraph %*% BN$marginal) > 0) * 1
		rownames(wholegraph) = rownames(BN$coregraph)
		wholegraph = wholegraph[colnames(BN$data),]
		for(reporter in 1:nrow(BN$data)){
			predicted.pos = sum(wholegraph[,reporter])
			if (predicted.pos == 0){predicted.pos = 1}
			parameters[3,reporter] = (BN$data[reporter,] %*% wholegraph[,reporter])/predicted.pos
			parameters[4,reporter] = sum(((BN$data[reporter,] - parameters[3,reporter])^2) %*% wholegraph[,reporter]) / predicted.pos
			predicted.neg = sum(1 - wholegraph[,reporter])
			if(predicted.neg == 0){predicted.neg = 1}
			parameters[1,reporter] = (BN$data[reporter,]%*%(1-wholegraph[,reporter]))/predicted.neg
			parameters[2,reporter] = sum(((BN$data[reporter,] - parameters[1,reporter])^2) %*% (1 - wholegraph[,reporter])) / predicted.neg
		}
		colnames(parameters) = rownames(BN$data)
		rownames(parameters) = c("mu0","var0","mu1","var1")
		return(parameters)
	}
}

parameters_continuous_Bayesian=function(BN, mu0 = 10, mu1 = 10, alpha1 = 10, beta1 = 1, alpha0 = 10, beta0 =1, v0 = 0.001, v1 = 0.001){
	parameters = matrix(1, nrow = 8, ncol = nrow(BN$data))
	wholegraph = ((BN$coregraph %*% BN$marginal) > 0) * 1
	rownames(wholegraph) = rownames(BN$coregraph)
	wholegraph = wholegraph[colnames(BN$data),]
	score = matrix(1, nrow = nrow(BN$data), ncol = 2)
	for(reporter in 1:nrow(BN$data)){
		N1 = sum(wholegraph[,reporter])
		div1 = N1
		if(div1 == 0){div1 = 1}
		mean1 = (BN$data[reporter,] %*% wholegraph[,reporter])/div1
		s1 = ((BN$data[reporter,] - mean1)^2) %*% wholegraph[,reporter]
		parameters[5, reporter] = alpha1 + N1
		parameters[6, reporter] = beta1 + s1 + v1*N1/(v1+N1) * (mean1 - mu1)^2
		parameters[7, reporter] = (v1*mu1 + N1*mean1)/(v1+N1)
		parameters[8, reporter] = v1 + N1
		N0 = sum(1 - wholegraph[,reporter])
		div0 = N0
		if(div0 == 0){div0 = 1}
		mean0 = (BN$data[reporter,]%*%(1-wholegraph[,reporter]))/div0
		s0 = ((BN$data[reporter,] - mean0)^2) %*% (1 - wholegraph[,reporter])
		parameters[1, reporter] = alpha0 + N0
		parameters[2, reporter] = beta0 + s0 + v0 * N0 / (v0 + N0) * (mean0 - mu0)^2
		parameters[3, reporter] = (v0 * mu0 + N0 * mean0)/(v0 + N0)
		parameters[4, reporter] = v0 + N0
	}
	colnames(parameters) = rownames(BN$data)
	rownames(parameters) = c("alpha0", "beta0", "mu0", "v0", "alpha1", "beta1", "mu1", "v1")
	return(parameters)
}


score_continuous_ML = function(BN, loc = NULL){
	if(!is.null(loc)){
		parameters = BN$parameters[,loc]
		wholegraph = ((BN$coregraph %*% BN$marginal[,loc]) > 0) * 1
		names(wholegraph) = rownames(BN$marginal)
		wholegraph = wholegraph[colnames(BN$data)]
		score1 = sum(dnorm(BN$data[loc, which(wholegraph == 1)], parameters[3], sqrt(parameters[4]), log = TRUE))
		score0 = sum(dnorm(BN$data[loc, which(wholegraph == 0)], parameters[1], sqrt(parameters[2]), log = TRUE))
		return(score1 + score0)
	}
	else{
		parameters = BN$parameters
		wholegraph = ((BN$coregraph %*% BN$marginal) > 0) * 1
		rownames(wholegraph) = rownames(BN$coregraph)
		wholegraph = wholegraph[colnames(BN$data),]
		score = matrix(1, nrow = nrow(BN$data), ncol = 2)
		for(reporter in 1:nrow(BN$data)){
			score[reporter, 1] = sum(dnorm(BN$data[reporter,which(wholegraph[,reporter] == 1)], parameters[3,reporter], sqrt(parameters[4,reporter]), log = TRUE))
			score[reporter, 2] = sum(dnorm(BN$data[reporter,which(wholegraph[,reporter] == 0)], parameters[1,reporter], sqrt(parameters[2,reporter]), log = TRUE))
		}
		return(sum(score))
	}
}

score_continuous_Bayesian = function(BN, loc = NULL, mu0 = 10, mu1 = 10, alpha1 = 10, beta1 = 1, alpha0 = 10, beta0 =1, v0 = 0.001, v1 = 0.001, likelihood = FALSE){
	if(!is.null(loc)){
		parameters = numeric(4)
		wholegraph = ((BN$coregraph %*% BN$marginal[,loc]) > 0) * 1
		names(wholegraph) = rownames(BN$marginal)
		wholegraph = wholegraph[colnames(BN$data)]
		N1 = sum(wholegraph)
		div1 = N1
		if(div1 == 0){div1 = 1}	
		parameters[3] = (BN$data[loc,] %*% wholegraph)/div1
		parameters[4] = sum(((BN$data[loc,] - parameters[3])^2) %*% wholegraph)
		N0 = sum(1 - wholegraph)
		div0 = N0
		if(div0 == 0){div0 = 1}
		parameters[1] = (BN$data[loc,] %*% (1 - wholegraph))/div0
		parameters[2] = ((BN$data[loc,] - parameters[1])^2) %*% (1 - wholegraph)
		score1 = (v1 / (v1 + N1)) ^ (1 / 2)  * gamma((alpha1 + N1) / 2) /  abs(beta1 + parameters[4] + v1 * N1 / (v1 + N1) * (parameters[3] - mu1) ^ 2) ^ ((alpha1 + mu1) / 2)
		score2 = (v0 / (v0 + N0)) ^ (1 / 2)  * gamma((alpha0 + N0) / 2) /  abs(beta0 + parameters[2] + v0 * N0 / (v0 + N0) * (parameters[1] - mu0) ^ 2) ^ ((alpha0 + mu0) / 2)
		score = score1 * score2
		return(log(score))
	}
	else{
		parameters = BN$parameters
		wholegraph = ((BN$coregraph %*% BN$marginal) > 0) * 1
		rownames(wholegraph) = rownames(BN$coregraph)
		wholegraph = wholegraph[colnames(BN$data),]
		score = matrix(1, ncol = 1, nrow = nrow(BN$data))
		for(reporter in 1:nrow(BN$data)){
			N1 = sum(wholegraph[,reporter])
			div1 = N1
			if(div1 == 0){div1 = 1}	
			parameters[3,reporter] = (BN$data[reporter,] %*% wholegraph[,reporter])/div1
			parameters[4,reporter] = sum(((BN$data[reporter,] - parameters[3,reporter])^2) %*% wholegraph[,reporter])
			N0 = sum(1 - wholegraph[,reporter])
			div0 = N0
			if(div0 == 0){div0 = 1}
			parameters[1,reporter] = (BN$data[reporter,]%*%(1-wholegraph[,reporter]))/div0
			parameters[2,reporter] = sum(((BN$data[reporter,] - parameters[1,reporter])^2) %*% (1 - wholegraph[,reporter]))
			score1 = (v1 / (v1 + N1)) ^ (1 / 2)  * gamma((alpha1 + N1)/2) /  abs(beta1 + parameters[4,reporter] + v1 * N1 / (v1 + N1) * (parameters[3,reporter] - mu1) ^ 2) ^ ((alpha1 + mu1) / 2)
			score2 = (v0 / (v0 + N0)) ^ (1 / 2)  * gamma((alpha0 + N0)/2) /  abs(beta0 + parameters[2,reporter] + v0 * N0 / (v0 + N0) * (parameters[1,reporter] - mu0) ^ 2) ^ ((alpha0 + mu0) / 2)
			score[reporter] = score1 * score2
		}
		if(likelihood){
			score = score * 1/pi^(ncol(BN$data) / 2) * abs(beta0) ^ (alpha0 / 2) * abs(beta1) ^ (alpha1 / 2) / (gamma(alpha1/2) * gamma(alpha0/2))
			}
		score = log(score)
		return(sum(score))
	}
}

# Input: a positive integer containing the number of bits
# Output: a sequence of positions which have to be altered successively to obtain a gray code
graychange = function(n)
{
    changepos = NULL
    for (j in 1:n)
    {
        changepos = c(changepos,j,changepos)
    }
    return(changepos)
}

exhaustive_BN = function (BN, verbose = TRUE, trace = FALSE, mode="binary_ML", IC.mode = "Bayesian", nu = c(0.001,0.1, 1,Inf), ...){	
	nrofedges = length(BN$variable.edges)
    changesequence = graychange(nrofedges)
    edgenr = BN$variable.edges[changesequence]
    original = BN$coregraph
    BN$coregraph = transitive.closure(BN$coregraph, mat=TRUE, loops=TRUE)
    diag(BN$coregraph) = diag(original)
    lastcoregraph = BN$coregraph
    BN$marginal = optimizemarginal(BN, mode, verbose = verbose, IC.mode = IC.mode, nu = nu, ...)
    if(mode == "binary_ML"){BN$parameters = parameters_discrete_ML(BN)
	    bestpost = score_discrete_ML(BN)}
    else if(mode== "binary_Bayesian"){bestpost = score_discrete_Bayesian(BN,...)}
	else if(mode== "continuous_Bayesian"){bestpost = score_continuous_Bayesian(BN,...)}
    else if(mode== "continuous_ML"){BN$parameters = parameters_continuous_ML(BN)
	    bestpost = score_continuous_ML(BN)}
	else{stop("inexisting mode")}
	BN$score = bestpost
	if(trace){post = rep(-Inf, 2^nrofedges)
		post[1] = bestpost
	}
    if(verbose){k = 1}
	if(verbose){cat("score: ", bestpost," ")}
	bestBN = BN
	for (j in 1:(2^nrofedges - 1)) {
		BN$coregraph = original
		BN$coregraph[edgenr[j]] = 1 - BN$coregraph[edgenr[j]]    
		if(is.dag(BN$coregraph)){
			original = BN$coregraph        
			BN$coregraph = transitive.closure(BN$coregraph, mat = TRUE, loops = TRUE)
    		diag(BN$coregraph) = diag(original)
    		if(identical(lastcoregraph, BN$coregraph)){next()}
    		if(verbose){k = k + 1}
    		lastcoregraph = BN$coregraph
    		BN$marginal = optimizemarginal(BN, mode, verbose = verbose, IC.mode = IC.mode, nu = nu, ...)
	        if(mode == "binary_ML"){
		        BN$parameters = parameters_discrete_ML(BN)
		        BN$score  = score_discrete_ML(BN)}
	    	else if(mode == "continuous_ML"){BN$parameters = parameters_continuous_ML(BN)
		    	BN$score  = score_continuous_ML(BN)}
		    else if(mode == "binary_Bayesian"){BN$score = score_discrete_Bayesian(BN,...)}
	    	else if(mode == "continuous_Bayesian"){BN$score = score_continuous_Bayesian(BN,...)}
	    	if(verbose){
	    		cat(" nr: ", k, "score: ", BN$score, " ")
	        	if(k %% 5 == 0){cat("\n")}}
	        if (BN$score > bestpost) {
	            bestpost = BN$score
	            bestBN   = BN
	        }
	        if(trace){post[j + 1] = BN$score}
		}
    }
    if(mode== "binary_Bayesian"){bestBN$parameters = parameters_discrete_Bayesian(BN,...)}
	else if(mode== "continuous_Bayesian"){bestBN$parameters = parameters_continuous_Bayesian(BN,...)}
	    	
    if (trace)
        return(list(bestBN,post))
    return(bestBN)
}

optimizemarginal = function(BN, mode, nu = c(0,0.001,1,10,Inf), trace = FALSE, verbose = FALSE, IC.mode = "Bayesian", ...){
	list = as.list(1:length(nu))
	IC = numeric(length(nu))
	nr=1
	for(penalty in nu){
		if(penalty == Inf){BN$marginal = optimizemarginal.nuInf(BN, mode, ...)}
		else{
			for (rep in 1:ncol(BN$marginal)){
				pr = numeric(4)
				BN$marginal[,rep]=0
				score.rep2 = -Inf
				continue = TRUE
				k = 1
				while(continue){
					score.rep1 = score.rep2
					poss.vertices = BN$poss.vertices[which(BN$marginal[BN$poss.vertices,rep] < 1)]
					if(length(poss.vertices) < 1){continue = FALSE
						break()}
					sc = rep(-Inf,times = nrow(BN$marginal))
					for (sign in poss.vertices){
						BN$marginal[sign,rep]=1	
						if(mode == "binary_ML"){
							BN$parameters[,rep] = parameters_discrete_ML(BN,rep)
							score = score_discrete_ML(BN,rep)
							}
						else if(mode == "continuous_ML"){
							BN$parameters[,rep] = parameters_continuous_ML(BN,rep)
							score = score_continuous_ML(BN,rep)
							}
						else if(mode == "binary_Bayesian"){score = score_discrete_Bayesian(BN,rep,...)}
						else if(mode == "continuous_Bayesian"){score = score_continuous_Bayesian(BN,rep,...)}
						else{stop("inexistant mode")}
						sc[sign] = score
						BN$marginal[sign,rep] = 0
					}
					score.rep2 = max(sc)
					new.edges = which(sc == max(sc))
					BN$marginal[new.edges,rep] = 1
					if(mode == "binary_ML"){
						BN$parameters[,rep] = parameters_discrete_ML(BN,rep)
						score.rep2 = score_discrete_ML(BN,rep)
						}
					else if(mode == "continuous_ML"){
						BN$parameters[,rep] = parameters_continuous_ML(BN,rep)
						score.rep2 = score_continuous_ML(BN,rep)
						}
					else if(mode == "binary_Bayesian"){score.rep2 = score_discrete_Bayesian(BN,rep,...)}
					else if(mode == "continuous_Bayesian"){score.rep2 = score_continuous_Bayesian(BN,rep,...)}
					else{stop("inexistant mode")}
				
					if(k == 1){continue = (score.rep2 - score.rep1 > (length(new.edges) - 1) * penalty)}
					else{continue = (score.rep2 - score.rep1 > length(new.edges) * penalty)}
					if(!continue){BN$marginal[new.edges,rep] = 0
						if(k == 1){BN$marginal[which.is.max(sc), rep] = 1}
					}
					k = k + 1
				}
			}
		}
		list[[nr]] = BN$marginal
		if(mode == "binary_ML"){
			BN$parameters = parameters_discrete_ML(BN)
			likelihood = score_discrete_ML(BN)
		}
		else if(mode== "binary_Bayesian"){likelihood = score_discrete_Bayesian(BN, likelihood = TRUE, ...)}
	    else if(mode== "continuous_Bayesian"){likelihood = score_continuous_Bayesian(BN, likelihood = TRUE, ...)}
		else if(mode== "continuous_ML"){BN$parameters = parameters_continuous_ML(BN)
	    	likelihood = score_continuous_ML(BN)}
	    if(IC.mode == "Bayesian"){	
			IC[nr] =  2 * (likelihood) - log(sum(BN$datatotal)) * sum(BN$marginal)}
		else if(IC.mode == "Akaike"){	
			IC[nr] =  2 * (likelihood) - 2 * sum(BN$marginal)}
		else (stop("inexistand IC.mode"))
			nr = nr + 1
	}
	if (verbose)
		{cat("nu.opt = ",nu[which.max(IC)]," nr marg not one connection = ", sum(apply(list[[which.max(IC)]],2,sum)!=1), " ")}
	if(trace){return(list(list,IC))}
	return(list[[which.max(IC)]])
}


optimizemarginal.nuInf = function(BN, mode, ...){
	Daten = BN$Daten
	for (rep in 1:ncol(BN$marginal)){
		sc = rep(-Inf,times = nrow(BN$marginal))
		pr = rep(1,times = 4)
		Daten = BN$Daten
		BN$marginal[,rep]=0
		for (sign in BN$poss.vertices){
			BN$marginal[sign,rep]=1	
				if(mode == "binary_ML"){
				BN$parameters[,rep] = parameters_discrete_ML(BN,rep)
				score = score_discrete_ML(BN,rep)
				}
			else if(mode == "continuous_ML"){
				BN$parameters[,rep] = parameters_continuous_ML(BN,rep)
				score = score_continuous_ML(BN,rep)
				}
			else if(mode == "binary_Bayesian"){score = score_discrete_Bayesian(BN,rep,...)}
			else if(mode == "continuous_Bayesian"){score = score_continuous_Bayesian(BN,rep,...)}
			else{stop("inexistant mode")}
			sc[sign] = score
			BN$marginal[sign,rep]=0
			}
		BN$marginal[which.is.max(sc),rep] = 1
		}
	graph = BN$marginal
	return(graph)
	}


ingreed_BN = function(BN, verbose = FALSE, mode = "binary_ML", trace = TRUE, max.iter = 100, IC.mode = "Bayesian", nu = c(0,0.001,1,0.1,Inf), ...) 
{
    if (trace) {
        tracelist = list()
    }
    if (verbose) 
        cat("Number of iterations : ")
    counter = 0
    repeat {
        coregraphold = BN$coregraph
        BN$marginal = optimizemarginal(BN, mode, verbose = verbose, IC.mode = IC.mode, nu = nu, ...)
        BN$coregraph = optimizecoregraph(BN, mode, ...)
        if (identical(BN$coregraph, coregraphold))
            {break()}
        counter = counter + 1
        if (verbose){
            cat("nr of iterations ", counter, " ")
            if ((counter%%5) == 0) 
                cat("\n")
        }
        if (trace)
            tracelist[[counter]] = BN
        if (counter == max.iter) 
            break()
    }
	if(mode == "binary_ML"){
        BN$parameters = parameters_discrete_ML(BN)
        BN$score = score_discrete_ML(BN)}
	else if(mode == "continuous_ML"){BN$parameters = parameters_continuous_ML(BN)
    	BN$score = score_continuous_ML(BN)}
    else if(mode == "binary_Bayesian"){BN$score = score_discrete_Bayesian(BN,...)
    	BN$parameters = parameters_discrete_Bayesian(BN,...)}
	else if(mode == "continuous_Bayesian"){BN$score = score_continuous_Bayesian(BN,...)
		BN$parameters = parameters_continuous_Bayesian(BN,...)}
	if (verbose) 
        cat("\n")
    if (!trace) 
        return(BN)
    return(list(BN = BN, tracelist = tracelist))
}

optimizecoregraph = function(BN, mode, best.marg=FALSE,...){
	bestcoregraph = BN$coregraph
	original = BN$coregraph
	original.red = transitive.reduction(original)
	if(mode == "binary_ML"){
		BN$parameters = parameters_discrete_ML(BN)
		bestpost = score_discrete_ML(BN)
		}
	else if(mode == "binary_Bayesian"){bestpost=score_discrete_Bayesian(BN,...)}
    else if(mode == "continuous_Bayesian"){bestpost=score_continuous_Bayesian(BN,...)}
	else if(mode == "continuous_ML"){BN$parameters = parameters_continuous_ML(BN)
	    bestpost = score_continuous_ML(BN)}
	else{stop("inexisting mode")}
    for (i in BN$variable.edges){
	    core = original.red
	    core[i] = 1 - core[i]
		BN$coregraph = transitive.closure(core, mat=TRUE, loops=TRUE)
		diag(BN$coregraph) = diag(original)
		if(!is.dag(core) || identical(original,BN$coregraph)){next()}
		if(best.marg){BN$marginal = optimizemarginal(BN, mode, verbose = verbose, IC.mode = IC.mode, nu = nu, ...)}
		if(mode == "binary_ML"){
			BN$parameters = parameters_discrete_ML(BN)
			post = score_discrete_ML(BN)
			}
		else if(mode == "binary_Bayesian"){post=score_discrete_Bayesian(BN,...)}
		else if(mode == "continuous_Bayesian"){post=score_continuous_Bayesian(BN,...)}
		else if(mode == "continuous_ML"){BN$parameters = parameters_continuous_ML(BN)
	    	post = score_continuous_ML(BN)}
	    else{stop("inexistant mode")}
    	if (post > bestpost){
			bestcoregraph = BN$coregraph
			bestpost=post}
		}
    return(bestcoregraph)
}

is.dag=function(graph){
	diag(graph) = 0
	graph1 = graph
	for (i in 1:nrow(graph))
		{graph1 = graph1 %*% graph
		if(sum(graph1)==0){return(TRUE)}
	}
	return(FALSE)
}

# simulate_BN = function(vert = 4, edges = NULL, gamma = 1, nr_intven = 3, beta = 0.9, reporters = 10, nullnode = TRUE, maxOutDegree = vert, maxInDegree = vert){
# 	core = diag(vert)
# 	colnames(core) = as.character(1:vert)
# 	rownames(core) = as.character(1:vert)	
# 	if(is.null(edges)){
# 		for(i in 1:vert){ # connect network randomly
# 			maxOutDegree = min( i,maxOutDegree)	
# 			degprob = (0:maxOutDegree)^(-gamma)
# 			degprob[1] = 1
# 			degprob = degprob/sum(degprob)
# 			outdeg = sample(0:maxOutDegree, 1, prob = degprob)# power law for out-degree => scale-free network
# 			if(outdeg > 0){
# 				idx0 = which(core[i,] == 0 & 1:vert < i)
# 				if(length(idx0) > 1){
# 					idx = which(colSums(core[,idx0]) <= maxInDegree)
# 					if(length(idx) > 0){			
# 						idx = sample(idx0[idx],min(outdeg,length(idx0[idx])),replace=TRUE)
# 						core[i,idx] = 1					
# 					}			
# 				}
# 				else if(length(idx0) == 1){
# 					idx = which(sum(core[idx0]) <= maxInDegree)
# 					if(length(idx) > 0){			
# 						idx = sample(idx0[idx],min(outdeg,length(idx0[idx])),replace=TRUE)
# 						core[i,idx] = 1					
# 					}			
# 				}	
# 			}
# 		}							
# 	}
# 	else{
# 		ve = which(core ==0)
# 		veind = which(core ==0, arr.ind = TRUE)
# 		veinddiff = veind[,1]-veind[,2]
# 		wh = which(veinddiff > 0)
# 		wh = ve[wh]
# 		sample = sample(wh,size = edges, replace = FALSE)
# 		core[sample]=1
# 	}
# 	original = core
#     core = transitive.closure(core, mat=TRUE, loops=TRUE)
#     diag(core) = diag(original)
# 	intven = rep(1:vert, each = nr_intven)	
# 	Delta = t(core)[,intven]
# 	ind = sample(1:nrow(core), reporters, replace = TRUE)
# 	Delta = Delta[ind,]
# 	colnames(Delta) = paste("I", intven,sep=".")
# 	if(nullnode){
# 		l = nrow(core)
# 		coregraph = matrix(0,nrow = l+1, ncol = l+1)
# 		coregraph[1:l,1:l] = core
# 		rownames(coregraph) = c(unique(colnames(Delta)),"null")
# 		colnames(coregraph) = c(unique(colnames(Delta)),"null")
# 		core = coregraph}
# 	Theta = matrix(0, ncol=reporters, nrow =ncol(core))
# 	for(i in 1:ncol(Theta))
# 		{Theta[ind[i],i] = 1}
# 	if(nullnode){rownames(Theta) = c(unique(colnames(Delta)),"null")}
# 	else{rownames(Theta) = unique(colnames(Delta))} 
# 	colnames(Theta) = as.character(1:ncol(Theta))
# 	report = rbinom(length(Delta), size = 1, prob = beta)
# 	report[which(Delta == 0)] = 1 - report[which(Delta == 0)]
# 	dim(report) = dim(Delta)
# 	colnames(report)=colnames(Delta)
# 	rownames(report) = colnames(Theta)
# 	BN = createBN(data = report, coregraph = core, marginal = Theta, nullnode = nullnode)
# 	return(BN)
# }
# 
# plotcoregraph = function (coregraph, main = NULL, transitive.reduction = TRUE) 
# {
#     gr = coregraph
#     diag(gr) = 0
#     wh1 = which(apply(gr,1,sum) > 0)
#     wh2 = which(apply(gr,2,sum) > 0)
#     wh = union(wh1,wh2)
#     if (length(wh)<1)(stop("coregraph empty"))
#     gr = gr[wh,wh]
#     if(transitive.reduction){gr = transitive.reduction(gr)}
#     grph = as(gr, "graphNEL")
#     plot(grph, main = main)
# }