
/*
matlab compilation options : 

 mex matlabcross.c -Isrc/include -Isrc/include/radiometer -Isrc/bin -lradiometer -Lbuild/lib  -L/home/michal/lscsoft/lal/lib -llal -lfftw3 -lfftw3f -lgsl -lgslcblas -lm -lgsl -lgslcblas -I/home/michal/lscsoft/lal/include/ 
*/


#include "mex.h"   
#include "stdio.h"
#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <getopt.h>
#include <fftw3.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_sf_coupling.h>
#include <radiometer/instrument.h>
#include <radiometer/sh_series.h>
#include <radiometer/inject.h>
#include <radiometer/correlator.h>
#include <instruments.h>
#include <stdlib.h>

static double wigner_3j(int ja, int jb, int jc, int ma, int mb, int mc)
{
	return gsl_sf_coupling_3j(2 * ja, 2 * jb, 2 * jc, 2 * ma, 2 * mb, 2 * mc);
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  
  double *conditionedData;
  double *windowData;
  int integrationLength, nTimeBins, segmentStartIndex ;
  double *outputArray;
  double *outI;
  double *Fp11;
  double *Fp12;
  double *Fp22;
  
  const int *dimArray;
  int i,k,iTimeBin,n,iCoeff;
  int minFreq,maxFreq;

  int p,kl,m;

  conditionedData=mxGetPr(prhs[0]);
  windowData=mxGetPr(prhs[1]);
  integrationLength = (int)mxGetScalar(prhs[2]);
  nTimeBins = (int)mxGetScalar(prhs[3]);
  segmentStartIndex = (int)mxGetScalar(prhs[4]);
  minFreq = (int)mxGetScalar(prhs[5]);
  maxFreq = (int)mxGetScalar(prhs[6]);
  Fp11=mxGetPr(prhs[7]);
  Fp12=mxGetPr(prhs[8]);
  Fp22=mxGetPr(prhs[9]);

  dimArray=mxGetDimensions(prhs[0]);
  int timeserieslength = dimArray[0];
  dimArray=mxGetDimensions(prhs[1]);
  double delta_t = 1.0 / 4096.0;

  /* define H1 and L1 : KLUDGE */
  double c = 299792458.0;
  /*    struct instrument *i1= instrument_new(4546374/c,
					842990/c,
					4378577/c);*/
  struct instrument *i1= instrument_new(-2.1614e6/c,
                                        -3.8347e6/c,
                                        4.6004e6/c);
  struct instrument *i2= instrument_new(-0.0743e6/c,
                                        -5.4963e6/c,
                                        3.2243e6/c);

  const struct instrument *instruments[] = {i1,i2};
  int n_instruments = sizeof(instruments) / sizeof(*instruments);
  
  /* compute spherical correlations plan */
  struct correlator_network_baselines *baselines = correlator_network_baselines_new(instruments, n_instruments);
  struct correlator_network_plan_td *tdplans = correlator_network_plan_td_new(baselines, delta_t);  
  struct correlator_network_plan_td *tdslowplans = correlator_network_plan_td_new(baselines, integrationLength*delta_t);  
  double *windows[] = {
    correlator_square_window_new(nTimeBins - 2 * tdplans->plans[0]->transient, 0, 1.0 / (nTimeBins - 2 * tdplans->plans[0]->transient))};
  
  
  
  const double *time_series_a=conditionedData;
  const double *time_series_b=conditionedData+timeserieslength;
  const double *window = windows[0];   
  struct correlator_plan_td *plan = tdplans->plans[0];
  struct correlator_plan_td *slowplan = tdslowplans->plans[0];
  
  /* allocate arrays for autocorrelation terms */ 
  
  int oneDlength = sh_series_length(plan->sample_a->l_max,1);
  int oneDslowlength = sh_series_length(slowplan->sample_a->l_max,1);
  int maxLength=oneDlength;
  plhs[1] = mxCreateDoubleScalar((double)maxLength);
  struct sh_series *autocor1 = sh_series_new_zero(plan->sample_a->l_max, 1);
  struct sh_series *autocor2 = sh_series_new_zero(plan->sample_a->l_max, 1);
  double *autot1 = malloc(integrationLength * sizeof(*autot1));
  double *autot2 = malloc(integrationLength * sizeof(*autot2));
  
  complex double *tmp_a=malloc(integrationLength * oneDlength * sizeof(*tmp_a));
  complex double *tmp_b=malloc(integrationLength * oneDlength * sizeof(*tmp_b));
  complex double *tmpauto1=malloc(integrationLength * oneDlength * sizeof(*tmpauto1));
  complex double *tmpauto2=malloc(integrationLength * oneDlength * sizeof(*tmpauto2));
  
  double *tmpTserie = malloc(nTimeBins * sizeof( *tmpTserie ));
  
  /* allocate output array, the array is filled with zeroes */

   plhs[0] = mxCreateDoubleMatrix(1, integrationLength*nTimeBins*maxLength,mxCOMPLEX); 
   /* plhs[0] = mxCreateDoubleMatrix(1, integrationLength*nTimeBins,mxREAL);  */
  outputArray = mxGetPr(plhs[0]);
  outI = mxGetPi(plhs[0]);

  /* precompute fft plan*/
  fftw_complex *in, *out;
  fftw_plan pp;
  in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * integrationLength);
  out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * integrationLength);
  pp = fftw_plan_dft_1d(integrationLength, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
  
  printf("%s\n","Pre-computing is finished");

  /* shift to the begining of the time series FIXME I may be wrong,
   create a sanity to check*/
  time_series_a=time_series_a + segmentStartIndex - 1 - plan->transient;
  time_series_b=time_series_b + segmentStartIndex - 1 - plan->transient;
  for(iTimeBin=0;iTimeBin < nTimeBins;iTimeBin++) { 
    
    for(k=0;k<integrationLength;k++){
      sh_series_array_dot(plan->sample_a, plan->proj_a, time_series_a++);
      sh_series_array_dot(plan->sample_b, plan->proj_b, time_series_b++);
      
      for(iCoeff=0;iCoeff<oneDlength;iCoeff++) {
	tmp_a[iCoeff*integrationLength + k] = plan->sample_a->coeff[iCoeff];
	tmp_b[iCoeff*integrationLength + k] = plan->sample_b->coeff[iCoeff];
	/*	tmpauto1[iCoeff*integrationLength + k] = 
	  autocor1->coeff[iCoeff];
	tmpauto2[iCoeff*integrationLength + k] = 
	autocor2->coeff[iCoeff];*/
      }
   
    }
    for(iCoeff=0;iCoeff<oneDlength;iCoeff++) {
      for(k=0;k<integrationLength;k++) {
	in[k] = tmp_a[iCoeff*integrationLength + k]*windowData[k];
      }
      fftw_execute(pp);
      for(k=0;k<integrationLength;k++) {
	tmp_a[iCoeff*integrationLength + k] = out[k];
	in[k] = tmp_b[iCoeff*integrationLength + k]*windowData[k];
      }
      fftw_execute(pp);
      for(k=0;k<integrationLength;k++) {
	tmp_b[iCoeff*integrationLength + k] = out[k];
      }  
    }
    /* ------------------------------------------------------- */
    for(k=minFreq-1;k<maxFreq;k++) {
   
      
      /*const int a_l_max = oneDlength;
      const int b_l_max = oneDlength;
      const int d_l_max = oneDlength;
      for(p = 0; p <= maxLength; p++) {

	/*	for(kl = 0; kl <oneDlength ; kl++) {
	  for(m = 0 ; m <oneDlength ; m++)
	    outputArray[iTimeBin*integrationLength + k +
			p*integrationLength*nTimeBins] += 
	      sqrt((2 * p + 1) * (2 * kl + 1) * (2 * m + 1) / (4 * M_PI)) * wigner_3j(kl, m, p, 0, 0, 0) * wigner_3j(kl, m, p, 0, 0, 0) *
	      creal(tmp_a[ kl*integrationLength + k])* 
	      ( creal(tmp_a[m*integrationLength + k]) /*- 
	      I*cimag(tmp_a[m*integrationLength + k])*//* );
							  }*//*
	outputArray[iTimeBin*integrationLength + k +
		    p*integrationLength*nTimeBins] += 
	  creal(tmp_a[ p*integrationLength + k]); 
	

	  }*/

      /*********/
      for(iCoeff=0;iCoeff<oneDlength;iCoeff++) {
	plan->sample_a->coeff[iCoeff] = 
	  tmp_a[iCoeff*integrationLength + k] ;
	plan->sample_b->coeff[iCoeff] = 
	  creal(tmp_a[iCoeff*integrationLength + k]) + 
	  -I*cimag(tmp_a[iCoeff*integrationLength + k]);
      }
      sh_series_product(plan->product, plan->sample_a, plan->sample_b, plan->product_plan);    
      for(iCoeff=0;iCoeff<oneDlength;iCoeff++) {
	outputArray[iTimeBin*integrationLength + k +
		    iCoeff*integrationLength*nTimeBins] += 
	  Fp11[k]*creal(plan->product->coeff[iCoeff]);
	outI[iTimeBin*integrationLength + k +
	     iCoeff*integrationLength*nTimeBins] += 
	  Fp11[k]*cimag(plan->product->coeff[iCoeff]);
	  }
      /**********/
      for(iCoeff=0;iCoeff<oneDlength;iCoeff++) {
	plan->sample_a->coeff[iCoeff] = 
	  tmp_a[iCoeff*integrationLength + k] ;
	plan->sample_b->coeff[iCoeff] = 
	  creal(tmp_b[iCoeff*integrationLength + k]) + 
	  -I*cimag(tmp_b[iCoeff*integrationLength + k]);
      }
      sh_series_product(plan->product, plan->sample_a, plan->sample_b, plan->product_plan);    
      for(iCoeff=0;iCoeff<oneDlength;iCoeff++) {
	outputArray[iTimeBin*integrationLength + k +
		    iCoeff*integrationLength*nTimeBins] += 
	  Fp12[k]*2*creal(plan->product->coeff[iCoeff]);
      }
      /*********/
        for(iCoeff=0;iCoeff<oneDlength;iCoeff++) {
	plan->sample_a->coeff[iCoeff] = 
	  tmp_b[iCoeff*integrationLength + k] ;
	plan->sample_b->coeff[iCoeff] = 
	  creal(tmp_b[iCoeff*integrationLength + k]) - 
	  I*cimag(tmp_b[iCoeff*integrationLength + k]) ;
      }
      sh_series_product(plan->product, plan->sample_a, plan->sample_b, plan->product_plan);    
      for(iCoeff=0;iCoeff<oneDlength;iCoeff++) {
	outputArray[iTimeBin*integrationLength + k +
		    iCoeff*integrationLength*nTimeBins] += 
	  Fp22[k]*creal(plan->product->coeff[iCoeff]);
	outI[iTimeBin*integrationLength + k +
	     iCoeff*integrationLength*nTimeBins] += 
	  Fp22[k]*cimag(plan->product->coeff[iCoeff]);
      }
      /* ----------------------------------------------------------- */
      /*for(iCoeff=0;iCoeff<oneDlength;iCoeff++) {
	plan->sample_a->coeff[iCoeff] = 
	  creal(tmp_b[iCoeff*integrationLength + k]) ;
	plan->sample_b->coeff[iCoeff] = 
	  creal(tmp_b[iCoeff*integrationLength + k]) ;
      }
      sh_series_product(plan->product, plan->sample_a, plan->sample_b, plan->product_plan);    
      for(iCoeff=0;iCoeff<oneDlength;iCoeff++) {
	outputArray[iTimeBin*integrationLength + k +
		    iCoeff*integrationLength*nTimeBins] += 
	  Fp22[k]*creal(plan->product->coeff[iCoeff]);
	outI[iTimeBin*integrationLength + k +
	     iCoeff*integrationLength*nTimeBins] += 
	  Fp22[k]*cimag(plan->product->coeff[iCoeff]);
	  }*/

      /*  for(iCoeff=0;iCoeff<oneDlength;iCoeff++) {
	plan->sample_a->coeff[iCoeff] = 
	  cimag(tmp_b[iCoeff*integrationLength + k]) ;
	plan->sample_b->coeff[iCoeff] = 
	  cimag(tmp_b[iCoeff*integrationLength + k]) ;
      }
      sh_series_product(plan->product, plan->sample_a, plan->sample_b, plan->product_plan);    
      for(iCoeff=0;iCoeff<oneDlength;iCoeff++) {
	outputArray[iTimeBin*integrationLength + k +
		    iCoeff*integrationLength*nTimeBins] += 
	  Fp22[k]*creal(plan->product->coeff[iCoeff]);
	outI[iTimeBin*integrationLength + k +
	     iCoeff*integrationLength*nTimeBins] += 
	  Fp22[k]*cimag(plan->product->coeff[iCoeff]);
	  }*/


      /* ------------------------------------------------------------ */

      /* for(iCoeff=0;iCoeff<oneDlength;iCoeff++) {
	
	plan->sample_a->coeff[iCoeff] = 
	  tmp_a[iCoeff*integrationLength + k] +
	  tmp_b[iCoeff*integrationLength + k] ;

	plan->sample_b->coeff[iCoeff] = 
	  creal(tmp_b[iCoeff*integrationLength + k]) + 
	  -I*cimag(tmp_b[iCoeff*integrationLength + k]) +
	  creal(tmp_a[iCoeff*integrationLength + k]) + 
	  -I*cimag(tmp_a[iCoeff*integrationLength + k]) ; 
	   
      }
      sh_series_product(plan->product, plan->sample_a, plan->sample_b, plan->product_plan);    
      for(iCoeff=0;iCoeff<oneDlength;iCoeff++) {
	outputArray[iTimeBin*integrationLength + k] += 
	  creal(plan->product->coeff[iCoeff]) * creal(plan->product->coeff[iCoeff]);
	  }*/
      
      

      /* ------------------------------------------------------------ */
      /* seperate sum squares */
      /* for(iCoeff=0;iCoeff<oneDlength;iCoeff++) {
	
	plan->sample_a->coeff[iCoeff] = 
	  tmp_b[iCoeff*integrationLength + k] ;

	plan->sample_b->coeff[iCoeff] = 
	  creal(tmp_a[iCoeff*integrationLength + k]) + 
	  -I*cimag(tmp_a[iCoeff*integrationLength + k]) ; 
	   
      }
      sh_series_product(plan->product, plan->sample_a, plan->sample_b, plan->product_plan);    
      for(iCoeff=0;iCoeff<oneDlength;iCoeff++) {
	outputArray[iTimeBin*integrationLength + k] += 
	  2*creal(plan->product->coeff[iCoeff])*creal(plan->product->coeff[iCoeff]) ;
      }
      for(iCoeff=0;iCoeff<oneDlength;iCoeff++) {
	plan->sample_a->coeff[iCoeff] = 
	  tmp_a[iCoeff*integrationLength + k] ;
	plan->sample_b->coeff[iCoeff] = 
	  creal(tmp_a[iCoeff*integrationLength + k]) -
	  I*cimag(tmp_a[iCoeff*integrationLength + k]) ; 
      }
      sh_series_product(plan->product, plan->sample_a, plan->sample_b, plan->product_plan);    
      for(iCoeff=0;iCoeff<oneDlength;iCoeff++) {
	outputArray[iTimeBin*integrationLength + k] += 
	  creal(plan->product->coeff[iCoeff])*creal(plan->product->coeff[iCoeff]);
      }
      
      for(iCoeff=0;iCoeff<oneDlength;iCoeff++) {
	plan->sample_a->coeff[iCoeff] = 
	  tmp_b[iCoeff*integrationLength + k] ;
	plan->sample_b->coeff[iCoeff] = 
	  creal(tmp_b[iCoeff*integrationLength + k]) - 
	  I*cimag(tmp_b[iCoeff*integrationLength + k]) ;
      }
      sh_series_product(plan->product, plan->sample_a, plan->sample_b, plan->product_plan);    
      for(iCoeff=0;iCoeff<oneDlength;iCoeff++) {
	outputArray[iTimeBin*integrationLength + k] += 
	  creal(plan->product->coeff[iCoeff])*creal(plan->product->coeff[iCoeff]);
      } 
      /* end seperate sum squares */

      /*    
      for(iCoeff=0;iCoeff<oneDlength;iCoeff++) {
	outputArray[iTimeBin*integrationLength + k +
		    iCoeff*integrationLength*nTimeBins] = 
	  creal(tmp_a[iCoeff*integrationLength + k]) ;
	outI[iTimeBin*integrationLength + k +
		    iCoeff*integrationLength*nTimeBins] = 
	  cimag(tmp_a[iCoeff*integrationLength + k]) ;
	
	  }*/
    
    }/* end loop over k (frequency) */
    
  }/* end loop ober iTimeBin */
  

  /*for(k=minFreq-1;k<maxFreq;k++) {
    for(iTimeBin=0;iTimeBin<nTimeBins;iTimeBin++) {
      tmpTserie[iTimeBin] = sqrt(outputArray[iTimeBin*integrationLength + k]);
    }
    for(iTimeBin=0;iTimeBin<nTimeBins-2*slowplan->transient;iTimeBin++) {
    sh_series_array_dot(slowplan->sample_a, slowplan->proj_a, tmpTserie+iTimeBin);

      outputArray[iTimeBin*integrationLength + k] = 0;
      for(iCoeff=0;iCoeff<oneDslowlength;iCoeff++) {
	outputArray[iTimeBin*integrationLength + k] +=
	  slowplan->sample_a->coeff[iCoeff]*slowplan->sample_a->coeff[iCoeff];
      }
    }
    } */
  
  correlator_network_plan_td_free(tdplans); 
  correlator_network_baselines_free(baselines);
  instrument_free(i2);  
  instrument_free(i1);  
  /* free(autot1);
     free(autot2);*/
  

  return;
}
    
