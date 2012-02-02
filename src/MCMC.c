/*
############################################################################
# Author:
# Date:
# Purpose:
# 
# 
# Related files: wrapper.c, simulation.R, MCMCdocu.txt, functionsMCMC.h
#		 netlearn.h, LogfilesMCMCrun.txt, and ouputMCMCdynoNEM
# 		 directory
#
#  
#      
############################################################################
*/
// Complete code to execute all the functions for the 
// metropolis hastings MCMC for dynoNEM
// the wrapper function for r is provided separately as wrapper.c file


#include "netlearn.h"


/*

// copy the initial network 
// modify the copied network
// check update factor if acceptable 
// calculate the MI with the older networks considering the timelags 
// count the number of accepted nework 
// if the number in greater than BURNIN and Mutual Information is less than acceptable
// continue modifying the last network
// calculate edge probability for all the networks in the sample

*/

//########## Function to calculate Hastings ratio /// 1. checked and changes suggested // ***** 2. Neighborhood removed

double updateFactor (double likLogOld, double logPriorOld, double logPriorScaleOld, double likLogNew, double logPriorNew, double logPriorScaleNew){ 
	return((likLogNew - likLogOld) + (logPriorNew - logPriorOld) + (logPriorScaleNew - logPriorScaleOld));   // hastings ratio ### Modified to handle Prior information		
}
// ###############################################
// Function to get log prior for the network 
//################################################
double logPrior(int nsgenes, double** net, double** prior, double inv_nu){
        int s;
        int k;

        if (prior==NULL) {
                return 0;
        }
        else {
                double p=0;
                for (s=0; s<nsgenes; s++) {
                        for (k=0; k<nsgenes; k++) {
                                p -= inv_nu*abs(net[s][k]-prior[s][k]);
                        }
                }
                p -= nsgenes*nsgenes*log(2);
                return p;
        }
        return 0;
}

// prior for 1/nu
double logPriorLambda(double inv_nu, double theta){
	return(log(theta) - theta * inv_nu);
}
//
//###########################################################################################################
//####################### log likelihood calculation (Plugged in from dynoNEM) #################### //Function updated

double** getPerturbProb(double** Psi, int T, int nsgenes, int k, double** perturb_prob){       
        int parent_perturb_prob;
        int s, t, p;        
        for(t = 0; t < T; t++){
                for(s = 0; s < nsgenes; s++){
                        perturb_prob[s][t] = 0; // default: s is unperturbed

                        perturb_prob[k][t] = 1; // the perturbed gene is always inactive

                        if(s != k){
                                for(p = 0; p < nsgenes; p++){
                                        if(Psi[p][s] != 0 && abs(Psi[p][s]) <= t){ // p is a parent
                                                if(t > 0){
                                                        parent_perturb_prob = perturb_prob[p][t-1];

                                                        }
                                                else
                                                        parent_perturb_prob = (p == k);
                                                if(parent_perturb_prob){
                                                        perturb_prob[s][t] = 1;
                                                        break;
                                                }
                                        }
                                }
                        }
                }
        
	}       	
	
	return(perturb_prob);
//
    		
//
}// working ###OK
// ########################################

double network_likelihood (double** Psi, int nsgenes, int negenes, int T, double*** D, double** egene_prior, int type, int nrep, double alpha, double beta, double*** perturb_prob){ // Egene_prior ????               
        int s, k, t, i;
    
        for(k = 0; k < nsgenes; k++){
                perturb_prob[k] = (double**) getPerturbProb(Psi, T, nsgenes, k, perturb_prob[k]);              
        }    
        	
        double tmp;
        double loglik0;
        double loglik = 0;
        
        double loglik_tmp;
        for (i=0; i<negenes; i++) {
                loglik_tmp=0;
                
                for (s=0; s < (nsgenes+1); s++) {
                        tmp=0;
                        for (k=0; k < nsgenes; k++) {
                                for (t=0; t < T; t++) {
                                        if(egene_prior[i][s] > 0 & s < nsgenes){
		                                  if(type == PVAL_DENS) // for p-value densities:
		                                        tmp += log((D[t][i][k]*perturb_prob[k][s][t] + (1-perturb_prob[k][s][t])*1) * egene_prior[i][s] + 1e-100); 
		                                  else if(type == EFFECT_PROB) // for effect probabilities
		                                        tmp += log((D[t][i][k]*perturb_prob[k][s][t] + (1-perturb_prob[k][s][t])*(1 - D[t][i][k])) * egene_prior[i][s] + 1e-100); 
		                                  else if(type == DISCRETE)// for count data:
		                                        tmp += log((pow(1-beta, D[t][i][k]*perturb_prob[k][s][t])*pow(beta, (nrep-D[t][i][k])*perturb_prob[k][s][t]) +
		                                               pow(alpha, D[t][i][k]*(1-perturb_prob[k][s][t]))*pow(1-alpha, (nrep-D[t][i][k])*(1-perturb_prob[k][s][t]))) * egene_prior[i][s]); 
					}
					if(egene_prior[i][s] > 0 & s == nsgenes){ // virtual "null" S-gene (not connected to any other S-gene, never perturbed) ==> automatic E-gene selection						
						if(type == PVAL_DENS)
							tmp += log(egene_prior[i][s] + 1e-100);
		                                else if(type == EFFECT_PROB) // for effect probabilities
		                                        tmp += log((1 - D[t][i][k]) * egene_prior[i][s] + 1e-100); 
		                                else if(type == DISCRETE) // for count data:
		                                        tmp += log((pow(alpha, D[t][i][k])*pow(1-alpha, (nrep-D[t][i][k]))) * egene_prior[i][s]);
					}	
					if(isnan(tmp) || !isfinite(tmp)){
						Rprintf("Numerical problem! tmp = NaN or Inf (i=%i, s=%i, k=%i, t=%i, D[t][i][k]=%g, egene_prior=%g)\n", i, s, k, t, D[t][i][k], egene_prior[i][s]);			
					}
                                }
                        }			
                        if (s==0) {
                                loglik0 = tmp;				
                        }
                        else {
                                loglik_tmp += exp(tmp - loglik0);
				if(isnan(loglik_tmp) || !isfinite(loglik_tmp)){
					Rprintf("Numerical problem! loglik_tmp = NaN or Inf (i=%i, s=%i, tmp=%g, loglik0=%g)\n", i, s, tmp, loglik0);			
				}
                        }
                }
                loglik += loglik0 +log(1 + loglik_tmp); // No O/P from here * removed b4 loglik     
        }
       return(loglik);
}

/*double** posteriorEGenePos(double** Psi, int nsgenes, int negenes, int T, double*** D, double** egene_prior, int type, int nrep, double alpha, double beta){
	double*** perturb_prob = (double***) calloc(nsgenes, sizeof(double**));
        int s, k, t, i;	
	for(i = 0; i < nsgenes; i++){
		perturb_prob[i] = (double**) calloc(nsgenes, sizeof(double*));// calloc changed to calloc
		for(k = 0; k < nsgenes; k++){
        	        perturb_prob[i][k] = (double*) calloc(T , sizeof(double));// calloc changed to calloc
			for(t = 0; t < T; t++)
				perturb_prob[i][k][t] = 0;
		}
        }
        for(k = 0; k < nsgenes; k++){
                perturb_prob[k] = (double**) getPerturbProb(Psi, T, nsgenes, k, perturb_prob[k]);
	}

        double tmp;      
	double** loglik_post = (double**) calloc(negenes, sizeof(double*));
        for (i=0; i<negenes; i++) {
		loglik_post[i] = (double*) calloc(nsgenes, sizeof(double));
                for (s=0; s<nsgenes; s++) {			
                        tmp=0.0;
                        for (k=0; k<nsgenes; k++) {
                                for (t=0; t<T; t++) {
                                        if(egene_prior[i][s] > 0){
                                          if(type == PVAL_DENS) // for p-value densities:
                                                tmp += log((D[t][i][k]*perturb_prob[k][s][t] + (1-perturb_prob[k][s][t])*1) * egene_prior[i][s]);
                                          else if(type == EFFECT_PROB) // for effect probabilities
                                                tmp += log((D[t][i][k]*perturb_prob[k][s][t] + (1-perturb_prob[k][s][t])*(1 - D[t][i][k])) * egene_prior[i][s]);
                                          else // for count data:
                                                tmp += log((pow(1-beta, D[t][i][k]*perturb_prob[k][s][t])*pow(beta, (nrep-D[t][i][k])*perturb_prob[k][s][t]) +
                                                       pow(alpha, D[t][i][k]*(1-perturb_prob[k][s][t]))*pow(1-alpha, (nrep-D[t][i][k])*(1-perturb_prob[k][s][t]))) * egene_prior[i][s]);
					}
                                }
                        }			
			loglik_post[i][s] = tmp;			
                }
        }
        for(k = 0; k < nsgenes; k++){
                for(s = 0; s < nsgenes; s++)
                        free(perturb_prob[k][s]);
                free(perturb_prob[k]);
        }	
        free(perturb_prob);
        return(loglik_post);
}*/

//############ Function to change the network // ***** To be reviewed // Modified
// Modified to remove the neighborhood ratio, random selection of altering operation
// Replaced the edge reversal operation with edge swpping operation

void alterNet(double** net, int nsgenes, int T, double** temp1){

  int i, j, noperations;
 
//   Rprintf("old network:\n");
//      for(i= 0; i<nsgenes; i++){
//     	for(j = 0; j<nsgenes; j++){
//     		Rprintf("%lf\t",net[i][j]);
//     	} Rprintf("\n");
//     }//
    
    copyNet(nsgenes, net, temp1); // ### OK
     
    //# 1) sample to indices i and j randomly
    //# 2) check which operations are possible
    //# 3) randomly choose any of the possible ones
    

    i = 0, j= 0;
    while (i == j){
      i = rand()%nsgenes;
      j = rand()%nsgenes;
      
    }
    //Rprintf("i = %d\tj = %d\n",i,j); //
    // # randomly pick index in array operations and execute the corresponding operation //
    // Generating random number
    
    int operations[3]; // changed to 4 to get edge deletion
    if (temp1[i][j] < T - 1) // increment operation
      operations[0] = 1;
    if(temp1[i][j] > 0) // decrement operation
      operations[1] = 1;
    if(temp1[i][j] != temp1[j][i]) // edge swap
      operations[2] = 1;
/*    if(temp1[i][j] > 0 && temp1[i][j] <= T) // delete operation
      operations[3] = 1;
 */   
    //

    noperations= 0;
    int k;
   
   for(k = 0; k < 3; k++){ // changed to 4 to get edge deletion else 3
      if(operations[k] == 1)
	noperations++;
    }  
   
    int r = rand()%noperations;
    //Rprintf("r = %d\n", r);

   //#$ Random number generated between 0 to noperations
    int kk = 0;
    for(k = 0; k < 3; k++){  // changed to 4 to get edge deletion else 3
      //#$ Rprintf("%d\n", operations[k]);
      if(operations[k] == 1){
// 	if(k == 0)
// 	  Rprintf("edge weight increment possible\n");
// 	 if(k == 1)
// 	  Rprintf("edge weight decrement possible\n");
// 	 if(k == 2)	   
// 	   Rprintf("edge weight swap possible\n");
	kk++;
      }
                     
      if(kk == r + 1){
	if(k == 0){ //# edge weight increment	
// 	 Rprintf("edge weight increment\n");
	 temp1[i][j] += 1;
	}
	else if(k == 1){ //# edge weight decrement
	  temp1[i][j] -= 1;
// 	  Rprintf("edge weight decrement\n");
	}    
	else if(k == 2){ //# edge weight swap
	  double tmp = temp1[i][j];
	  temp1[i][j] = temp1[j][i];
	  temp1[j][i] = tmp;
// 	  Rprintf("edge weight swap\n");
	}
// 	else if(k == 3){ //# edge delete // implemented again
// 	  temp1[i][j] = 0;
// 	}	// till this
	break;
      }	
    }
    
//       Rprintf("new network:\n");
//      for(i= 0; i<nsgenes; i++){
//     	for(j = 0; j<nsgenes; j++){
//     		Rprintf("%lf\t",temp1[i][j]);
//     	} Rprintf("\n");
//     }//

}  // ### OK (returns a new nwtwork) // 
//################## Mutual Information Functions
//
/*double mutInfo(double** net1, double** net2, int nsgenes, int T){

	double Pa, Pb, Pab, MI;

	int i, j, u, v;
	MI = 0.0;
	for(i = 0; i < T; i++){
		u = counter(net1, i, nsgenes);
		Pa = (double) u/ (double)(nsgenes*nsgenes);
		for(j = 0; j< T; j++){
			v = counter(net2, j, nsgenes);
			Pb = (double)v/(double)(nsgenes*nsgenes);
			Pab = jointProb(net1, net2, nsgenes, i, j);
			if (Pab != 0.0 && (Pa*Pb != 0.0))
				MI += Pab*log(Pab/(Pa*Pb));				
		}
	}	
	return (MI);
}

int counter(double** v1, int x, int nsgenes){
	int i, j;
	int count;
	count = 0;
	for(i=0; i<nsgenes; i++){
		for(j=0; j<nsgenes; j++){
			if(v1[i][j]==x)	
			count++;
		}
	}
return(count);
}

double jointProb(double** v1, double** v2, int nsgenes, int y, int x){
	int i, j; 
	int jointcount;
	jointcount = 0;    
	for (i = 0; i < nsgenes; i++){
		for(j = 0; j < nsgenes;j++){
			if((v2[i][j]==y) && (v1[i][j]==x)){
				jointcount++;
			}
		}
	}
	return ((double)jointcount / (double)(nsgenes*nsgenes));
}
*/
// ### OK calculates the mutual Information between two networks
//#######################################################################
//############//Function to copy network  

void copyNet(int nsgenes, double** net, double** netCopy)
{
	//printf("Copying Network\n");
	int i, j;
	for (i = 0; i < nsgenes; i++){
		for (j = 0; j < nsgenes; j++){
			netCopy[i][j] = net[i][j];
		}
	}
}
//############################ ###################

void MCMCrun(long sample, long burnin, double** net, int nsgenes, int negenes, int T, double*** D, double** networkPrior,  double **Egene_prior, double priorScale, double theta, int type, int nrep, double alpha, double beta, int seed, double* allLikelihoods, double** sdMat, double** matrix_r){ // Egene_prior ????
  Rprintf("SAMPLE = %ld\nBURNIN = %ld\nNSGENES = %d\nNEGENES = %d\nT = %d\nTYPE = %d\nNREP = %d\nALPHA = %lf\nBETA = %lf\nTHETA = %lf\n", sample, burnin, nsgenes, negenes, T, type, nrep, alpha, beta, theta); // ###OK
  //
    int i, j, t;
   //# // variable for likelihood sum initialization
    double loglikMean, loglikSum, mutinf, delta, logPrior_cur_scale;  
    srand(seed);
   //# // Sum of likelihoods 
    double** oldNet = (double**) calloc(nsgenes, sizeof(double*)); 
    double** matrix = (double**) calloc(nsgenes, sizeof(double*)); 	// for DIC calculation
//    double** matrix_r = (double**) calloc(nsgenes, sizeof(double*));	// rounded off matrix
    double** M = (double**) calloc(nsgenes, sizeof(double*)); 		// for variance calculation
    //double** varMat= (double**) calloc(nsgenes, sizeof(double*)); 	// for variance calculation
    double** var_mean = (double**) calloc(nsgenes, sizeof(double*));
    double** newNet = (double**) calloc(nsgenes, sizeof(double*));
   /* double*** networks = (double***) calloc(TIMELAG, sizeof(double**));           
    for(i = 0; i < TIMELAG; i++){ // BURNIN changed to nsgenes
	networks[i] = (double**) calloc(nsgenes, sizeof(double*));
	for(j = 0; j < nsgenes; j++) // SAMPLE changed to nsgenes
	  networks[i][j] = (double*) calloc (nsgenes, sizeof(double));
    }*/
    double*** perturb_prob = (double***) calloc(nsgenes, sizeof(double**));        
    for(i = 0; i < nsgenes; i++){ // changed to nsgenes from SAMPLE
	matrix[i] = (double*) calloc (nsgenes, sizeof(double));	// for DIC calculation
	//matrix_r[i] = (double*) calloc (nsgenes, sizeof(double));	// rounded off matrix
	M[i] = (double*)calloc(nsgenes, sizeof(double));		// for variance calculation
	//varMat[i]= (double*)calloc(nsgenes, sizeof(double));		// for variance calculation
	var_mean[i] = (double*) calloc(nsgenes, sizeof(double));
	newNet[i] = (double*) calloc (nsgenes, sizeof(double));	
	oldNet[i] = (double*) calloc (nsgenes, sizeof(double));	
	
	perturb_prob[i] = (double**) calloc(nsgenes, sizeof(double*));// calloc changed to calloc
	for(j = 0; j < nsgenes; j++){
		matrix[i][j] = 0;
		matrix_r[i][j] = 0;
		M[i][j] = 0;
		var_mean[i][j] = 0;
		newNet[i][j] = 0;
		oldNet[i][j] = 0;
                perturb_prob[i][j] = (double*) calloc(T , sizeof(double));// calloc changed to calloc
		for(t = 0; t < T; t++)
			perturb_prob[i][j][t] = 0;
        }
	
    }   
    // M
    long counter = 0;
    long stored = 0;
    long stored2 = 0;
    long nsampled = 0;
    int converged = 0;
    double hfactor;
    double likelihood, logPrior_cur, delta_loglambda, priorScale_new;    
    //
    Rprintf("counter = %ld and converged = %d \n",counter,converged);// ### OK
    //
    copyNet(nsgenes, net, oldNet);
    
    double likLogOld = network_likelihood(oldNet, nsgenes, negenes, T, D, Egene_prior, type, nrep, alpha, beta, perturb_prob);
    double logPriorOld = logPrior(nsgenes, oldNet, networkPrior, priorScale);
    double logPriorScale = logPriorLambda(priorScale, theta);
    long accept = 0;
    loglikSum = 0.0;     
    allLikelihoods[0] = likLogOld;
    GetRNGstate();
    //
   // copyNet(nsgenes, oldNet, networks[0]);
    while (counter < (burnin + sample)){	
	/*if(stored2 >= (TIMELAG-1) && !converged && counter < burnin){	     // stored2 >= TIMELAG    
	    mutinf = mutInfo(networks[TIMELAG-1], networks[0], nsgenes, T);  // networks[TIMELAG-1]
	    printf("The mutual information = %lf\n", mutinf);
	    // Update Convergence   
	    converged = (mutinf < THRESHOLD);	
	    if(converged){ // Stop burnin if converged and go to sampling
	      printf("Converged at %ld and Mutual Information is %lf\n", counter, mutinf);
	      counter=burnin;
	    }
	   
 	}*/
	if(counter % 100 == 0){
		delta_loglambda = rnorm(0, sqrt(0.5)); // sample such that new log-lambda is with 95% probability within one log-unit apart from current lambda
		priorScale_new *= pow(2, delta_loglambda);
		priorScale_new += 1e-7;
		logPrior_cur_scale = logPriorLambda(priorScale_new, theta);
	}
	else{
		priorScale_new = priorScale;		
		logPrior_cur_scale = logPriorScale;
	}
	alterNet(net, nsgenes, T, newNet);
	//
	likelihood = network_likelihood(newNet, nsgenes, negenes, T, D, Egene_prior, type, nrep, alpha, beta,  perturb_prob);	
	logPrior_cur = logPrior(nsgenes, newNet, networkPrior, priorScale_new);
		
	hfactor = updateFactor(likLogOld, logPriorOld, logPriorScale, likelihood, logPrior_cur, logPrior_cur_scale);
		
	//#// random number between 0 and 1 
	double r2 = 0;
	while(r2 == 0)
		r2 = ((double)(rand() % 100000001)) /((double)100000000);
	// Random number between 0 and 1 generated to be compared to the Hastings ratio
	if(hfactor >= log(r2)){	   // if hfactor is > 0, this condition is ALWAYS fullfilled	 (-INFINITY < log(r2) <= 0)    	    
	    copyNet(nsgenes, newNet, net);
	    priorScale = priorScale_new;
	    if(counter % 100 == 0)
		    Rprintf("new prior scale = %g\n", priorScale);
	    accept++;
	    
	  /*  if(!converged && counter <= burnin){ 		
		if(stored2 < TIMELAG - 1)
		    stored2++;
		else{ // shifts the array frame for networks one block backwards
		    for(i = 0; i < TIMELAG - 1; i++)
			copyNet(nsgenes, networks[i + 1], networks[i]);
		}		  
		copyNet(nsgenes, net, networks[stored2]);
	    }*/
      
	    likLogOld = likelihood;
	    logPriorOld = logPrior_cur;
	    logPriorScale = logPrior_cur_scale;
	}
        allLikelihoods[counter + 1] = likelihood; // Put likelihoods into an array
	if(counter % 100 == 0){
	  Rprintf("iter = %ld, accepted = %ld, likelihood = %g\n",counter, accept, likLogOld);	   // Number of Accepted networks 
	}
	counter++;	
					     	
    
        if((counter > burnin) && counter%100 == 0){
	//printf("Burnin (convergence) complete %ld and Mutual Information is %lf\n", burnin, mutinf);
		loglikSum += likelihood;

		nsampled++;
		//#// Storing likelihood sum
		for(i= 0; i < nsgenes; i++){
			    for(j = 0; j < nsgenes; j++){
				  matrix[i][j] += net[i][j]; // matrix updated 
				  delta = net[i][j]-var_mean[i][j];
				  var_mean[i][j] += delta/nsampled;
				  M[i][j] += (delta*(net[i][j]-var_mean[i][j]));
			    }
		 }
       }
   }
    PutRNGstate();
    Rprintf("\n\nnsampled = %ld\n", nsampled);
    Rprintf("Likelihood sum is %lf\n",loglikSum);
    loglikMean = loglikSum/((double)nsampled);
    Rprintf("Mean Likelihood is %lf\n",loglikMean);  
    Rprintf("SDs for the edges in network\n");
    for(i = 0; i < nsgenes; i++){      
	    for(j = 0; j < nsgenes; j++){
		matrix_r[i][j] = round(matrix[i][j]/((double)nsampled)) ;		// averaged and then rounded off ('round up' use ceil(), 'round down' use floor())
	    	sdMat[i][j] = sqrt(M[i][j]/((double)(nsampled - 1)));						// Variance calculation from M-matrix
		Rprintf("%lf\t",sdMat[i][j]);
	    }
	    Rprintf("\n");
    }      
      
      //Dhat;
   //  Rprintf("Converged at %ld and Mutual Information is %lf \n ", burnin, mutinf);
     double Dhat = network_likelihood(matrix_r, nsgenes, nsgenes, T, D, Egene_prior, type, nrep, alpha, beta, perturb_prob); // Egene_prior ????
     Rprintf("The Dhat is %lf\n",Dhat);
     double DIC = Dhat - 2*loglikMean;
     Rprintf("DIC is %lf\n", DIC);

     for(i = 0; i < nsgenes; i++){
	free(matrix[i]);
	free(M[i]);
	free(var_mean[i]);
	free(newNet[i]);
	free(oldNet[i]);
	for(j = 0; j < nsgenes; j++)
		free(perturb_prob[i][j]);
	free(perturb_prob[i]);
     }
     free(matrix);
     free(M);
     free(var_mean);
     free(newNet);
     free(oldNet);
     free(perturb_prob);
}
