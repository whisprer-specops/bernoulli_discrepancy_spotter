/* fastlapreg.h */
/* definitions for fast lapreg period search */
/*
[Adapated from.  Los Alamos National Security, LLC.]
 
 
 This program is free software; you can redistribute it 
 and/or modify it under the terms of the GNU General Public License 
 as published by the Free Software Foundation; version 2.0 of the License. 
 Accordingly, this program is distributed in the hope that it will be useful, 
 but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY 
 or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License 
 for more details
*/


extern char* version_fastlaplacereg;

#define FFTSIGNCONVENTION -1
#define MAXORDER    15

#endif


#ifndef _FASTLAPLACEREG_H_
#define _FASTLAPLACEREG_H_

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multimin.h>

  //  #include <gsl/gsl_errno.h>

#include <memory.h>
#include <string.h>

extern char* version_fastchi2;

#define FFTSIGNCONVENTION -1	/* In GSL, a forward FFT has exp( - 2 pi n/N) */
								/* Numerical Recipes would have   +  so +1    */

#define MAXORDER    15          /* both the detrending order and the harmonic order */
                                /* Can be set to any value */

/* Find the lapreg value for deviation from the mean */
/* if pmean == NULL, then deviation from 0 */

double laplacereg_objective(const gsl_vector *v, void *params);
void laplacereg_gradient(const gsl_vector *v, void *params, gsl_vector *grad);
void laplacereg_fdf(const gsl_vector *v, void *params, double *f, gsl_vector *grad);

/* Read the next source in the infile, and return the number of */
/* measurements (0 for err or eof) */
int     readsource(FILE* infile, 
                char *echo, double *sampletimes, float*vals, float *stderrs);

			/* Produce a vector of lapreg values as a function of frequency */
			/* dov: _d_ata _o_ver _v_ariance (direct or FFT of) */
			/* oov: _o_ne  _o_ver _v_ariance (direct or FFT of) */
void    fastLap(const double *sampletimes, const float *vals, const float *stderrs, 
                long ndata, long nharmonics, double deltat, 
                long nrealpoints, int nIterations, double tstart,
                float *dovstorage/*[nrealpoints]*/, float *oovstorage/*[nrealpoints]*/,
                float *lapreductionlist);
    
/* Calculate the lap-reg at a single frequency, returning the harmonic */
/* coefficients as well */
double	singleLap(const double *sampletimes, const float *vals, const float *stderrs, 
                long ndata, long nharmonics, double tstart,
                double frequency, float *harmoniccoeffs);
                
/* Optimize the lap-reg near a single frequency +/- fhwidth,   */
/* returns the chisquare improvement, the harmonic coeffs, and the  */
/* optimal frequency */
double	localOptimizeLap(const double *sampletimes, const float *vals, const float *stderrs, 
                  long ndata, long nharmonics, double tstart,
                  double frequency, double fhwidth, 
                  float *harmoniccoeffs, double *bestf);


/* Calculate the height of the parabolic maximum near p0, if any, */
/* and zero if not peak.  This is not a general routine: */
/* it is only useful if the values are all positive and you don't */
/* care what the answer is except near the maximum */
inline float parabolic_array_peak(float* p0);

/* Return the index of the best lap-reg in the array, where 'best' may */
/* be adjusted by the number of trials. */
int    findBestGridPoint(float *lapreductionvalues, 
                         int findexmin, int findexmax,
                         int nharmonics,
                         int fAdjusted);


double findBestFrequency(int fAdjusted,
                         const double *sampletimes, const float *vals, const float *stderrs, 
                         long ndata, int nharmonics,
                         double freqmin, double deltaf, double t0, float lapmargin,
                         float *lapreductionvalues, int nlapvalues, 
                         float *pbest_peakup);

/* given the set of FFT values for data/variance and 1/variance     */
/* calculate the best lapreg reduction.  The gsl_{matrix,vector}s must be */
/* pre-allocated for efficiency, and are used to solve: */
/*      alpha beta = x      */
/* cos*[0] is the mean value (f=0),	sin*[0] is ignored, */
/* *data[0..nharmonics], *var[0..2*nharmonics]  */
double	calcLapReduction(int nharmonics, 
                const float *cosdata, const float *sindata,
                const float *cosvar, const float *sinvar,
				gsl_matrix *alpha, gsl_vector *beta, gsl_vector *x);
                
            /* Calculate the alpha and beta matrices */
void  calcAlphaBeta(int nharmonics, 
                const float *cosdata, const float *sindata,
                const float *cosvar, const float *sinvar,
                gsl_matrix *alpha, gsl_vector *beta);


double  Mfunc(int m_index, double phi);


/* Return a value which is monotonic in the probability of achieving */
/* the given lapreg reduction by chance, adjusted for the number */
/* of trials if desired.  This allows relative comparison of two laplace */
/* probabilities (which must have the same ndof and fAdjusted) without having */
/* to worry about e.g. Prob(lapreg) being zero due to underflow. */
double Pmetric(double lapreg_reduction, double ndof, int ntrials, int fAdjusted);

/* return the chiqsuared reduction required for this pmetric value */
double Pmetric_inverse(double p, double ndof, int ntrials, int fAdjusted);


/* utility functions */
double dmin(double d1, double d2);
double dmax(double d1, double d2);


#endif  /*  _FASTLAPLACEREG_H_ */
