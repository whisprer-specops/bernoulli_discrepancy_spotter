/* fastchi2.h */
/* definitions for fast chi2 period search */

#ifndef _FASTCHI2_H_
#define _FASTCHI2_H_

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_errno.h>
#include <memory.h>
#include <string.h>

extern char* version_fastchi2;

#define FFTSIGNCONVENTION -1	/* In GSL, a forward FFT has exp( - 2 pi n/N) */
								/* Numerical Recipes would have   +  so +1    */

#define MAXORDER    15          /* both the detrending order and the harmonic order */
                                /* Can be set to any value */

/* Find the chisquared value for deviation from the mean */
/* if pmean == NULL, then deviation from 0 */
double chisqconst(int ndata, float *values, float *stderrs, float* pmean);

/* Read the next source in the infile, and return the number of */
/* measurements (0 for err or eof) */
int     readsource(FILE* infile, 
                char *echo, double *sampletimes, float*vals, float *stderrs);

			/* Produce a vector of chisquared values as a function of frequency */
			/* dov: _d_ata _o_ver _v_ariance (direct or FFT of) */
			/* oov: _o_ne  _o_ver _v_ariance (direct or FFT of) */
void    fastChi(const double *sampletimes, const float *vals, const float *stderrs, 
                long ndata, long nharmonics, double deltat, 
                long nrealpoints, long nchivalues, double tstart,
                float *dovstorage/*[nrealpoints]*/, float *oovstorage/*[nrealpoints]*/,
                float *chireductionlist);
    
/* Calculate the chi-squared at a single frequency, returning the harmonic */
/* coefficients as well */
double	singleChi(const double *sampletimes, const float *vals, const float *stderrs, 
                long ndata, long nharmonics, double tstart,
                double frequency, float *harmoniccoeffs);
                
/* Optimize the chi-squared near a single frequency +/- fhwidth,   */
/* returns the chisquare improvement, the harmonic coeffs, and the  */
/* optimal frequency */
double	localOptimizeChi(const double *sampletimes, const float *vals, const float *stderrs, 
                  long ndata, long nharmonics, double tstart,
                  double frequency, double fhwidth, 
                  float *harmoniccoeffs, double *bestf);


/* Calculate the height of the parabolic maximum near p0, if any, */
/* and zero if not peak.  This is not a general routine: */
/* it is only useful if the values are all positive and you don't */
/* care what the answer is except near the maximum */
inline float parabolic_array_peak(float* p0);

/* Return the index of the best chi-squared in the array, where 'best' may */
/* be adjusted by the number of trials. */
int    findBestGridPoint(float *chireductionvalues, 
                         int findexmin, int findexmax,
                         int nharmonics,
                         int fAdjusted);


double findBestFrequency(int fAdjusted,
                         const double *sampletimes, const float *vals, const float *stderrs, 
                         long ndata, int nharmonics,
                         double freqmin, double deltaf, double t0, float chimargin,
                         float *chireductionvalues, int nchivalues, 
                         float *pbest_peakup);

/* given the set of FFT values for data/variance and 1/variance     */
/* calculate the best chisquared reduction.  The gsl_{matrix,vector}s must be */
/* pre-allocated for efficiency, and are used to solve: */
/*      alpha beta = x      */
/* cos*[0] is the mean value (f=0),	sin*[0] is ignored, */
/* *data[0..nharmonics], *var[0..2*nharmonics]  */
double	calcChiReduction(int nharmonics, 
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
/* the given chi-squared reduction by chance, adjusted for the number */
/* of trials if desired.  This allows relative comparison of two chisquared */
/* probabilities (which must have the same ndof and fAdjusted) without having */
/* to worry about e.g. Prob(chisquared) being zero due to underflow. */
double Pmetric(double chisq_reduction, double ndof, int ntrials, int fAdjusted);

/* return the chiqsuared reduction required for this pmetric value */
double Pmetric_inverse(double p, double ndof, int ntrials, int fAdjusted);


/* utility functions */
double dmin(double d1, double d2);
double dmax(double d1, double d2);

#endif  /*  _FASTCHI2_H_ */
