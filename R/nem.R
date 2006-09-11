nem <- function(D,inference="pairwise",models=NULL,type="mLL",para=NULL,hyperpara=NULL,Pe=NULL,Pm=NULL,local.prior.size=length(unique(colnames(D))),local.prior.bias=1,verbose=TRUE){

#------------------------------
# Sanity checks                

if (!(inference %in% c("pairwise","search")))   stop("\nnem> argument 'inference' is not valid\n")
if (!(type %in% c("mLL","FULLmLL")))            stop("\nnem> argument 'type' is not valid")
if (is.null(para)   & is.null(hyperpara))       stop("\nnem> provide either 'para' or 'hyperpara'\n")
if (type=="mLL"     & is.null(para))            stop("\nnem> provide argument 'para'\n")
if (type=="FULLmLL" & is.null(hyperpara))       stop("\nnem> provide argument 'hyperpara'\n")
if (!is.null(hyperpara)){
    if (length(hyperpara)!=4)                   stop("\nnem> 'hyperpara' is not a vector of length 4")
    if (!all(hyperpara > 0))                    stop("\nnem> 'hyperpara' must be >0")
    }
if (!is.null(para)){          
    if (length(para)!=2)                        stop("\nnem> 'para' is not a vector of length 2")
    if (any(para < 0) | any(para > 1))          stop("\nnem> 'para' must be in [0,1]")
    }


Sgenes <- unique(colnames(D))


#------------------------------
# PAIRWISE                     

if (inference == "pairwise"){
        if (local.prior.size <= 0 | local.prior.bias <= 0) stop("\nnem> local prior parameters invalid")
        Pm <- local.model.prior(local.prior.size,length(Sgenes),local.prior.bias)
        result <- pairwise.posterior(D,type,para,hyperpara,Pe,Pm,verbose)     
    }


#------------------------------
# SEARCH                       

if (inference == "search"){ 
        if (is.null(models)) models <- enumerate.models(length(Sgenes),Sgenes,verbose)
        result <- score(models,D,type,para,hyperpara,Pe,verbose)
       }


#------------------------------
# OUTPUT                       
#class(result) <- "nem"
return(result)

}
