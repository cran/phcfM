////////////////////////////////////////////////////////////////////
// demography_fixed.cc
//
// demography.cc samples from the posterior distribution of a
// Gaussian linear regression model
//
////////////////////////////////////////////////////////////////////
//
// Original code by Ghislain Vieilledent, March 2012
// CIRAD UR BSEF
// ghislain.vieilledent@cirad.fr / ghislainv@gmail.com
//
////////////////////////////////////////////////////////////////////
//
// This software is distributed under the terms of the GNU GENERAL
// PUBLIC LICENSE Version 2, June 1991.  See the package LICENSE
// file for more information.
//
// Copyright (C) 2011 Ghislain Vieilledent
// 
////////////////////////////////////////////////////////////////////
//
// Revisions: 
// - G. Vieilledent, on March 2012
//
////////////////////////////////////////////////////////////////////


// Scythe libraries
#include "matrix.h"
#include "distributions.h"
#include "stat.h"
#include "la.h"
#include "ide.h"
#include "smath.h"
#include "rng.h"
#include "mersenne.h"
// R libraries
#include <R.h> // needed to use Rprintf()
#include <R_ext/Utils.h> // needed to allow user interrupts

using namespace scythe;
using namespace std;

extern "C"{

    /* Gibbs sampler function */
    void demography_fixed (
	
	// Constants and data
	const int *ngibbs, const int *nthin, const int *nburn, // Number of iterations, burning and samples
	const int *nobs, // Constants
	const int *np, // Number of fixed and random covariates
	const double *Y_vect, // Observed response variable
	const double *X_vect, // Covariate for fixed effects
        // Parameters to save
	double *beta_vect, // Fixed effects
	double *V, // Variance of residuals
	// Defining priors
	const double *mubeta_vect, const double *Vbeta_vect,
	const double *s1_V, const double *s2_V,
	// Diagnostic
	double *Deviance,
	double *Y_pred, // Fitted values (predictive posterior mean)
	// Seeds
	const int *seed,
	// Verbose
	const int *verbose
	
	) {
	
	////////////////////////////////////////////////////////////////////////////////
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	// Defining and initializing objects
	
        ///////////////////////////
	// Redefining constants //
	const int NGIBBS=ngibbs[0];
	const int NTHIN=nthin[0];
	const int NBURN=nburn[0];
	const int NSAMP=(NGIBBS-NBURN)/NTHIN;
	const int NOBS=nobs[0];
	const int NP=np[0];

        ///////////////
	// Constants //

	// Small fixed matrices indexed on i for data access
	Matrix<double> *Xi_arr = new Matrix<double>[NOBS];
	Matrix<double> *Yi_arr = new Matrix<double>[NOBS];
	for(int i=0; i<NOBS; i++) {
	    Xi_arr[i] = Matrix<double>(1,NP);
	    Yi_arr[i] = Matrix<double>(1,1);
	    for (int p=0; p<NP; p++) {
		Xi_arr[i](p)=X_vect[p*NOBS+i];
	    }
	    Yi_arr[i](0)=Y_vect[i];
	    Y_pred[i]=0; // We initialize Y_pred to zero to compute the predictive posterior mean
	} 

	////////////////////
	// Priors objects //
	Matrix<double> mubeta(NP,1,mubeta_vect);
	Matrix<double> Vbeta(NP,NP,Vbeta_vect);

	/////////////////////////////////////
	// Initializing running parameters //
	Matrix<double> beta_run(NP,1,false); // Unicolumn matrix of fixed effects
	for (int p=0; p<NP; p++) {
	    beta_run(p)=beta_vect[p*NSAMP];
	}
	double V_run=V[0];
	double Deviance_run=Deviance[0];

	////////////
	// Message//
	Rprintf("\nRunning the Gibbs sampler. It may be long, keep cool :)\n\n");
	R_FlushConsole();
	//R_ProcessEvents(); for windows	

	/////////////////////
	// Set random seed //
	mersenne myrng;
	myrng.initialize(*seed);

	///////////////////////////////////////////////////////////////////////////////////////
	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	// Gibbs sampler

	for (int g=0;g<NGIBBS;g++) {

	    //////////////////////////////////
	    // vector beta: Gibbs algorithm //

	    // invVi, sum_V and sum_v
	    Matrix<double> sum_V(NP,NP);
	    Matrix<double> sum_v(NP,1);
	    for (int i=0; i<NOBS; i++) {
	    	sum_V += crossprod(Xi_arr[i]);
	    	sum_v += t(Xi_arr[i])*Yi_arr[i];
	    }

	    // big_V
	    Matrix<double> big_V=invpd(invpd(Vbeta)+sum_V/V_run);
	    
	    // small_v
	    Matrix<double> small_v=invpd(Vbeta)*mubeta+sum_v/V_run;

	    // Draw in the posterior distribution
	    beta_run=myrng.rmvnorm(big_V*small_v,big_V);


	    ////////////////
	    // variance V //
	    
	    // e
	    Matrix<double> e(1,1,true,0.0);
	    for (int i=0; i<NOBS; i++) {
	    	e+=crossprod(Yi_arr[i]-Xi_arr[i]*beta_run);
	    }

	    // Parameters
	    double S1=*s1_V+(NOBS/2); //shape
	    double S2=*s2_V+0.5*e(0); //rate

            // Draw in the posterior distribution
	    V_run=1/myrng.rgamma(S1,S2);

	    //////////////////////////////////////////////////
	    //// Deviance

	    // logLikelihood
	    double logLk=0;
	    for (int i=0; i<NOBS; i++) {
	    	// Y_hat
		Matrix<double> Y_hat=Xi_arr[i]*beta_run;
	    	// L
		logLk+=log(dnorm(Y_vect[i],Y_hat(0),sqrt(V_run)));
	    }

	    // Deviance
	    Deviance_run=-2*logLk;


	    //////////////////////////////////////////////////
	    // Output
	    if(((g+1)>NBURN) && (((g+1)%(NTHIN))==0)) {
	    	int isamp=((g+1)-NBURN)/(NTHIN);
		for (int p=0; p<NP; p++) {
		    beta_vect[p*NSAMP+(isamp-1)]=beta_run(p);
		}
		V[isamp-1]=V_run;
		Deviance[isamp-1]=Deviance_run;
		for (int i=0; i<NOBS; i++) {
		    Matrix<double> Y_hat=Xi_arr[i]*beta_run;
		    Y_pred[i]+=Y_hat(0)/NSAMP; // We compute the predictive posterior mean of NSAMP values 
		}
	    }
	    
	    
	    //////////////////////////////////////////////////
	    // Progress bar
	    double Perc=100*(g+1)/(NGIBBS);
	    if(((g+1)%(NGIBBS/100))==0 && (*verbose==1)){  
	    	Rprintf("*");
	    	R_FlushConsole();
	    	//R_ProcessEvents(); for windows
	    	if(((g+1)%(NGIBBS/10))==0){
		    Rprintf(":%.1f%%\n",Perc);
      	    	    R_FlushConsole();
	    	    //R_ProcessEvents(); for windows
	    	}
	    } 


            //////////////////////////////////////////////////
	    // User interrupt
	    R_CheckUserInterrupt(); // allow user interrupts 
	    
	} // Gibbs sampler

	///////////////
	// Delete memory allocation (see new)
	delete[] Yi_arr;
	delete[] Xi_arr;

    } //  end demography_fixed function

} // end extern "C"

////////////////////////////////////////////////////////////////////
// END
////////////////////////////////////////////////////////////////////
