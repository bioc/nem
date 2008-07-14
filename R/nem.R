nem <- function(D,inference="nem.greedy",models=NULL,type="mLL",para=NULL,hyperpara=NULL,Pe=NULL,Pm=NULL,Pmlocal=NULL,local.prior.size=length(unique(colnames(D))),local.prior.bias=1,triples.thrsh=0.5,lambda=0,delta=1,selEGenes=FALSE,trans.close=TRUE,verbose=TRUE){
#------------------------------
# Sanity checks              
if(inference == "nem.greedyMAP")
	type = "CONTmLLMAP"
  
if (!(inference %in% c("pairwise", "triples", "search","ModuleNetwork","RelNet","nem.greedy","nem.greedyMAP"))) 
	stop("\nnem> argument 'inference' is not valid\n")
if (!(type %in% c("mLL", "FULLmLL", "CONTmLL","CONTmLLBayes","CONTmLLMAP", "CONTmLLDens", "CONTmLLRatio"))) 
	stop("\nnem> argument 'type' is not valid")
if (is.null(para) & is.null(hyperpara) & !(type %in% c("CONTmLL", "CONTmLLBayes", "CONTmLLMAP", "CONTmLLDens", "CONTmLLRatio")))
	stop("\nnem> provide either 'para' or 'hyperpara'\n")
if (inference != "RelNet" & type == "mLL" & is.null(para)) 
	stop("\nnem> provide argument 'para'\n")
if (type == "FULLmLL" & is.null(hyperpara)) 
	stop("\nnem> provide argument 'hyperpara'\n")
if (!is.null(hyperpara)) {
	if (length(hyperpara) != 4) 
		stop("\nnem> 'hyperpara' is not a vector of length 4")
	if (!all(hyperpara > 0)) 
		stop("\nnem> 'hyperpara' must be >0")
}
if (!is.null(para)) {
	if (length(para) != 2) 
		stop("\nnem> 'para' is not a vector of length 2")
	if (any(para < 0) | any(para > 1)) 
		stop("\nnem> 'para' must be in [0,1]")
}
if(lambda < 0) lambda <- abs(lambda)

Sgenes <- unique(colnames(D))
if(selEGenes){	
	return(nem.featureselection(D, inference, models, type, para, hyperpara, Pe, Pm, Pmlocal, local.prior.size, local.prior.bias, triples.thrsh, lambda, trans.close=trans.close, verbose))
}

#------------------------------
# PAIRWISE                     

if (inference == "pairwise"){	
	if (is.null(local.prior.size))
		Pmlocal <- NULL 
	else{ 
		if (local.prior.size <= 0 | local.prior.bias <= 0) 
			stop("\nnem> local prior parameters invalid")
        	Pmlocal <- local.model.prior(local.prior.size,length(Sgenes),local.prior.bias)		
        }
        result <- pairwise.posterior(D,type,para,hyperpara,Pe,Pmlocal,Pm,lambda,delta,verbose)    
    }

#------------------------------
# MODULE NETWORK                       
else if(inference == "ModuleNetwork"){
	result <- moduleNetwork(D,type,Pe,Pm,lambda,delta,para,hyperpara,verbose=verbose)	
}

#------------------------------
# TRIPLES                     

else if (inference == "triples"){

        result <- triples.posterior(D,type,para,hyperpara,Pe,Pmlocal,Pm,lambda,delta,triples.thrsh,verbose)
        #A.t <- transitive.lp(A)

    }

#------------------------------
# GREEDY                     
else if(inference == "nem.greedy"){
	result <- nem.greedy(D,initial=models,type=type,Pe=Pe,Pm=Pm,lambda=lambda,delta=delta,para=para,hyperpara=hyperpara, verbose=verbose)	
}

else if(inference == "nem.greedyMAP"){
	print(delta)
	
	result <- nem.greedyMAP(D,Pe=Pe,Pm=Pm,lambda=lambda, delta=delta, trans.close=trans.close, verbose=verbose)	
}

#------------------------------
# SEARCH                       

else if (inference == "search"){ 
        if (is.null(models)) models <- enumerate.models(length(Sgenes),Sgenes,verbose)
        result <- score(models,D,type,para,hyperpara,Pe,Pm,lambda,delta,verbose)
}
else
	stop(paste("Unknown inference method", inference,"\n"))

#------------------------------
# OUTPUT                       
#class(result) <- "nem"
return(result)

}
