////////////////////////////////////////////////////////////////////
// deforestation.cc
//
// deforestation.cc samples from the posterior distribution of a
// logistic regression model with variable time-interval between censuses
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
#include <R.h>           // needed to use Rprintf()
#include <R_ext/Utils.h> // needed to allow user interrupts
// My own libraries
#include "logit_invlogit.h"

using namespace scythe;
using namespace std;

/////////////////////////////////////////////////////////////
// log(LikelihoodxPrior) function for Metropolis algorithm

static double logLk(const Matrix<double>& Y, const Matrix<double>& X, const Matrix<double>& Int, 
		    const Matrix<double>& beta) {
    // Likelihood
    Matrix<double> logit_theta = X*beta;
    Matrix<double> theta = 1.0 / (1.0 + exp(-logit_theta));
    double log_Lk = 0;
    for (int i=0; i<X.rows(); i++) {
	double theta_prim = 1-pow(1-theta(i),Int(i));
	log_Lk +=  Y(i)*log(theta_prim) + (1-Y(i))*log(1-theta_prim);
    }
    
    // Result
    return log_Lk;
}

/////////////////////////////////////////////////////////////
// log(Prior) function for Metropolis algorithm

static double logP(const Matrix<double>& beta, const Matrix<double>& mubeta, const Matrix<double>& Vbeta) {
    // Prior
    double log_P = 0.0;
    if (det(Vbeta)>0.0) {
	log_P = lndmvn(beta, mubeta, Vbeta);
    }
  
    // Result
    return log_P;
}

/////////////////////////////////////////////////////////////
// C++ function                        

extern "C"{
  
    /* Gibbs sampler function */
    void deforestation(

	// Constants and data
	const int *ngibbs, const int *nthin, const int *nburn, // Number of iterations, burning and samples
	const int *nobs, // Constants
	const int *np, // Number of fixed effects
	const int *Y_vect, // Observed response variable
	const double *X_vect, // Covariates
	const double *Int_vect, // Time-interval
	// Objects for proposal in Metropolis
	double *tune_scalar, // Tuning scalar: will change with adaptive algorithm
	const double *VCV_vect, // Variance-Covariance matrix for fixed effects from glm
        // Parameters to save
	double *beta_vect, // Fixed effects
	// Defining priors
	const double *mubeta_vect, const double *Vbeta_vect,
	// Diagnostic
	double *Deviance,
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

	// Data //
	Matrix<double> X (NOBS,NP);
	for(int i=0; i<NOBS; i++) {
	    for (int p=0; p<NP; p++) {
		X(i,p)=X_vect[p*NOBS+i];
	    }
	} 
	Matrix<double> Y (NOBS,1,Y_vect);
	Matrix<double> Int (NOBS,1,Int_vect);

	// Priors objects //
	const Matrix<double> mubeta(NP,1,mubeta_vect);
	const Matrix<double> Vbeta(NP,NP,Vbeta_vect);
	//Rprintf("%.1f",mubeta(0));
	//Rprintf("%.1f",Vbeta(0));

        // Objects for proposal in Metropolis //
	double tune_scalar_run = tune_scalar[0];
	Matrix<double> Id_NP = eye<double>(NP);
	const Matrix<double> VCV(NP,NP,VCV_vect);
	const Matrix<double> Vbeta_with_VCV = invpd(invpd(Vbeta)+invpd(VCV));

	/////////////////////////////////////
	// Initializing running parameters //
	Matrix<double> beta_run(NP,1,false); // Unicolumn matrix of fixed effects
	for (int p=0; p<NP; p++) {
	    beta_run(p)=beta_vect[p*NSAMP];
	}
	double Deviance_run=Deviance[0];

        ///////////////////////////////////////////
	// Acceptance rate for adaptive sampling //
	int nA=0; // Number of acceptance
	double Ar=0; // Acceptance rate

	////////////
	// Message//
	Rprintf("\nRunning the Gibbs sampler. It may be long, please keep cool :)\n\n");
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


            //////////////////////////////////////////////////
	    //// beta

	    // logLP_now
	    double logLk_now=logLk(Y, X, Int, beta_run);
	    double logLP_now=logLk_now+logP(beta_run, mubeta, Vbeta);

            // logLP_prop
	    Matrix<double> tune_run=tune_scalar_run*Id_NP;
	    Matrix<double> propV = tune_run * Vbeta_with_VCV * tune_run;
	    Matrix<double> propC = cholesky(propV);
	    Matrix<double> beta_prop=gaxpy(propC, myrng.rnorm(NP, 1, 0, 1), beta_run);
	    double logLP_prop=logLk(Y, X, Int, beta_prop)+logP(beta_prop, mubeta, Vbeta);

	    // ratio
	    double r=exp(logLP_prop-logLP_now); // ratio
	    double z=myrng.runif();

	    // Actualization
	    if (z < r) {
	    	beta_run=beta_prop;
	    	nA++;
	    }


	    //////////////////////////////////////////////////
	    //// Deviance
	    Deviance_run=-2*logLk_now;


	    //////////////////////////////////////////////////
	    // Output
	    if(((g+1)>NBURN) && (((g+1)%(NTHIN))==0)){
		int isamp=((g+1)-NBURN)/(NTHIN);
		for (int p=0; p<NP; p++) {
		    beta_vect[p*NSAMP+(isamp-1)]=beta_run(p);
		}
		Deviance[isamp-1]=Deviance_run;
	    }


	    ///////////////////////////////////////////////////////
	    // Adaptive sampling (on the burnin period)
	    int DIV=0;
	    if (NGIBBS >=1000) DIV=100;
	    else DIV=NGIBBS/10;
	    if((g+1)%DIV==0 && (g+1)<=NBURN){
		const double ropt=0.44;
		Ar=(double)nA/(double)DIV;
		if(Ar>=ropt) tune_scalar_run=tune_scalar_run*(2-(1-Ar)/(1-ropt));
		else tune_scalar_run=tune_scalar_run/(2-Ar/ropt);
		nA=0.0; // We reinitialize the number of acceptance to zero
	    }
	    if((g+1)%DIV==0 && (g+1)>NBURN){
		Ar=(double)nA/(double)DIV;
		nA=0.0; // We reinitialize the number of acceptance to zero
	    }

    
	    //////////////////////////////////////////////////
	    // Progress bar
	    double Perc=100*(g+1)/(NGIBBS);
	    if(((g+1)%(NGIBBS/100))==0 && verbose[0]==1) {  
	    	Rprintf("*");
	    	R_FlushConsole();
	    	//R_ProcessEvents(); for windows
	    	if(((g+1)%(NGIBBS/10))==0){
		    Rprintf(":%.1f%%, mean accept. rate=%.3f\n",Perc,Ar);
      	    	    R_FlushConsole();
	    	    //R_ProcessEvents(); for windows
	    	}
	    } 


            //////////////////////////////////////////////////
	    // User interrupt
	    R_CheckUserInterrupt(); // allow user interrupt 	    
	
	} // Gibbs sampler

	///////////////
	// Output for tune_scalar
	tune_scalar[0]=tune_scalar_run;
	
    } // end entry_mortality_gibbs_fixed function

} // end extern "C"

////////////////////////////////////////////////////////////////////
// END
////////////////////////////////////////////////////////////////////
