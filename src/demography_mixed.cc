////////////////////////////////////////////////////////////////////
// demography_mixed.cc
//
// demography_mixed.cc samples from the posterior distribution of a
// Gaussian hierarchical linear regression model
//
// The code uses Algorithm 2 of Chib & Carlin (1999) for efficient
// inference of (\beta | Y, sigma^2, Vb).
//
// Chib, S. & Carlin, B. P. (1999) On MCMC sampling in hierarchical
// longitudinal models, Statistics and Computing, 9, 17-26
//
////////////////////////////////////////////////////////////////////
//
// Original code by Ghislain Vieilledent, March 2012
// CIRAD UR B&SEF
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
    void demography_mixed (
	
	// Constants and data
	const int *ngibbs, const int *nthin, const int *nburn, // Number of iterations, burning and samples
	const int *nobs, const int *ngroup, // Constants
	const int *np, const int *nq, // Number of fixed and random covariates
	const int *IdentGroup, // Vector of group
	const double *Y_vect, // Observed response variable
	const double *X_vect, // Covariate for fixed effects
	const double *W_vect, // Covariate for random effects
        // Parameters to save
	double *beta_vect, // Fixed effects
	double *b_vect, // Random effects
	double *Vb_vect, // Variance of random effects
	double *V, // Variance of residuals
	// Defining priors
	const double *mubeta_vect, const double *Vbeta_vect,
	const double *r, const double *R_vect,
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
	const int NGROUP=ngroup[0];
	const int NP=np[0];
	const int NQ=nq[0];

	///////////////
	// Constants //
	// Number of observations by group k
	int *nobsk = new int[NGROUP];
	for (int k=0; k<NGROUP; k++) {
	    nobsk[k]=0;
	    for (int n=0; n<NOBS; n++) {
		if (IdentGroup[n]==k) {
		    nobsk[k]+=1;
		}
	    }
	}
	// Position of each group in the data-set
	int **posk_arr = new int*[NGROUP];
	for (int k=0; k<NGROUP; k++) {
	    posk_arr[k] = new int[nobsk[k]];
	    int repk=0;
	    for (int n=0; n<NOBS; n++) {
	    	if (IdentGroup[n]==k) {
	    	    posk_arr[k][repk]=n;
	    	    repk++;
	    	}
	    }
	}
	// Small fixed matrices indexed on k for data access
	Matrix<double> *Yk_arr = new Matrix<double>[NGROUP];
	Matrix<double> *Xk_arr = new Matrix<double>[NGROUP];
        Matrix<double> *Wk_arr = new Matrix<double>[NGROUP];
	for(int k=0; k<NGROUP; k++) {
	    Xk_arr[k] = Matrix<double>(nobsk[k],NP);
	    Wk_arr[k] = Matrix<double>(nobsk[k],NQ);
	    Yk_arr[k] = Matrix<double>(nobsk[k],1);
	    for (int m=0; m<nobsk[k]; m++) {
		for (int p=0; p<NP; p++) {
		    Xk_arr[k](m,p)=X_vect[p*NOBS+posk_arr[k][m]];
		}
		for (int q=0; q<NQ; q++) {
		    Wk_arr[k](m,q)=W_vect[q*NOBS+posk_arr[k][m]];
		}
		Yk_arr[k](m,0)=Y_vect[posk_arr[k][m]];
		Y_pred[posk_arr[k][m]]=0; // We initialize Y_pred to zero to compute the predictive posterior mean
	    }
	} 
	Matrix<double> *tXk_arr = new Matrix<double>[NGROUP];
	Matrix<double> *tWk_arr = new Matrix<double>[NGROUP];
	Matrix<double> *cpXk_arr = new Matrix<double>[NGROUP];
	Matrix<double> *tXWk_arr = new Matrix<double>[NGROUP];
	Matrix<double> *tWXk_arr = new Matrix<double>[NGROUP];
	Matrix<double> *tXYk_arr = new Matrix<double>[NGROUP];
	Matrix<double> *tWYk_arr = new Matrix<double>[NGROUP];
	for(int k=0; k<NGROUP; k++) {
	    tXk_arr[k] = t(Xk_arr[k]);
	    tWk_arr[k] = t(Wk_arr[k]);
	    cpXk_arr[k] = crossprod(Xk_arr[k]);
	    tXWk_arr[k] = t(Xk_arr[k])*Wk_arr[k];
	    tWXk_arr[k] = t(Wk_arr[k])*Xk_arr[k];
	    tXYk_arr[k] = t(Xk_arr[k])*Yk_arr[k];
	    tWYk_arr[k] = t(Wk_arr[k])*Yk_arr[k];
	}

	////////////////////
	// Priors objects //
	Matrix<double> mubeta(NP,1,mubeta_vect);
	Matrix<double> Vbeta(NP,NP,Vbeta_vect);
	Matrix<double> R(NQ,NQ,R_vect);

	/////////////////////////////////////
	// Initializing running parameters //
	Matrix<double> *bk_run = new Matrix<double>[NGROUP]; // Random effects
	for (int k=0;k<NGROUP;k++) {
	    bk_run[k] = Matrix<double>(NQ,1);
	    for (int q=0; q<NQ; q++) { 
		bk_run[k](q)=b_vect[q*NGROUP*NSAMP+k*NSAMP];
	    }
	}
	Matrix<double> beta_run(NP,1,false); // Unicolumn matrix of fixed effects
	for (int p=0; p<NP; p++) {
	    beta_run(p)=beta_vect[p*NSAMP];
	}
	Matrix<double> Vb_run(NQ,NQ,true,0.0);
	for (int q=0; q<NQ; q++) {
	    for (int qprim=0; qprim<NQ; qprim++) {
		Vb_run(q,qprim)=Vb_vect[qprim*NQ*NSAMP+q*NSAMP];
	    }
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
	// Gibbs sampler (see Chib et Carlin, 1999 p.4 "Blocking for Gaussian mixed models")

	for (int g=0;g<NGIBBS;g++) {

	    //////////////////////////////////
	    // vector beta: Gibbs algorithm //

	    // invVi, sum_V and sum_v
	    // see http://en.wikipedia.org/wiki/Woodbury_matrix_identity with A=(1/V_run)*Id
	    Matrix<double> sum_V(NP,NP);
	    Matrix<double> sum_v(NP,1);
	    for (int k=0; k<NGROUP; k++) {
	    	sum_V += (1/V_run)*cpXk_arr[k]-pow(1/V_run,2)*tXWk_arr[k]*invpd(invpd(Vb_run)+tWk_arr[k]*(1/V_run)*Wk_arr[k])*tWXk_arr[k];
	    	sum_v += (1/V_run)*tXYk_arr[k]-pow(1/V_run,2)*tXWk_arr[k]*invpd(invpd(Vb_run)+tWk_arr[k]*(1/V_run)*Wk_arr[k])*tWYk_arr[k];
	    }

	    // big_V
	    Matrix<double> big_V=invpd(invpd(Vbeta)+sum_V/V_run);
	    
	    // small_v
	    Matrix<double> small_v=invpd(Vbeta)*mubeta+sum_v/V_run;

	    // Draw in the posterior distribution
	    beta_run=myrng.rmvnorm(big_V*small_v,big_V);


            ///////////////////////////////
	    // vector b: Gibbs algorithm //

	    // Loop on group
	    for (int k=0; k<NGROUP; k++) {
    
	    	// big_Vk
	    	Matrix<double> big_Vk=invpd(invpd(Vb_run)+crossprod(Wk_arr[k])/V_run);

	    	// small_vk
	    	Matrix<double> small_vk=(t(Wk_arr[k])*(Yk_arr[k]-Xk_arr[k]*beta_run))/V_run;

	    	// Draw in the posterior distribution
	    	bk_run[k]=myrng.rmvnorm(big_Vk*small_vk,big_Vk);
	    }

	 
	    ////////////////////////////////////////////
	    // vector of variance Vb: Gibbs algorithm //
	    Matrix<double> SSBb(NQ,NQ,true,0.0);      
	    for(int k=0; k<NGROUP; k++) {
	    	SSBb+=bk_run[k]*t(bk_run[k]); 
	    }      
	    int Vb_dof=(*r)+(NGROUP);
	    Matrix<double> Vb_scale = invpd(SSBb+(*r)*R);   
	    Vb_run=invpd(myrng.rwish(Vb_dof,Vb_scale)); 


	    ////////////////
	    // variance V //
	    
	    // e
	    Matrix<double> e(1,1,true,0.0);
	    for (int k=0; k<NGROUP; k++) {
	    	e+=crossprod(Yk_arr[k]-Xk_arr[k]*beta_run-Wk_arr[k]*bk_run[k]);
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
	    for (int k=0; k<NGROUP; k++) {
	    	for (int m=0; m<nobsk[k]; m++) {
	    	    // Which one ?
	    	    int w=posk_arr[k][m];
	    	    // Y_hat
	    	    Matrix<double> Y_hat=Xk_arr[k](m,_)*beta_run+Wk_arr[k](m,_)*bk_run[k];
	    	    // L
	    	    logLk+=log(dnorm(Y_vect[w],Y_hat(0),sqrt(V_run)));
	    	}
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
		for (int k=0; k<NGROUP; k++) {
		    for (int q=0; q<NQ; q++) {
			b_vect[q*NGROUP*NSAMP+k*NSAMP+(isamp-1)]=bk_run[k](q);
		    }
		}
		for (int q=0; q<NQ; q++) {
		    for (int qprim=0; qprim<NQ; qprim++) {
			Vb_vect[qprim*NQ*NSAMP+q*NSAMP+(isamp-1)]=Vb_run(q,qprim);
		    }
		}
		V[isamp-1]=V_run;
		Deviance[isamp-1]=Deviance_run;
		for (int k=0; k<NGROUP; k++) {
		    for (int m=0; m<nobsk[k]; m++) {
			// Which one ?
			int w=posk_arr[k][m];
			// Y_hat
			Matrix<double> Y_hat=Xk_arr[k](m,_)*beta_run+Wk_arr[k](m,_)*bk_run[k];
			// Y_pred
			Y_pred[w]+=Y_hat(0)/NSAMP;
		    }
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
	delete[] nobsk;
	for(int k=0; k<NGROUP; k++) {
	    delete[] posk_arr[k];
	}
	delete[] posk_arr;
	delete[] Yk_arr;
	delete[] Xk_arr;
	delete[] Wk_arr;
	delete[] tXk_arr;
	delete[] tWk_arr;
	delete[] cpXk_arr;
	delete[] tXWk_arr;
	delete[] tWXk_arr;
	delete[] tXYk_arr;
	delete[] tWYk_arr;
	delete[] bk_run;

    } //  end demography_mixed function

} // end extern "C"

////////////////////////////////////////////////////////////////////
// END
////////////////////////////////////////////////////////////////////

