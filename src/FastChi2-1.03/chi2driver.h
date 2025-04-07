/* chi2driver.h */
/* definitions for driver of fastchi2 program */

/*
 Copyright 2007.  Los Alamos National Security, LLC. This material was 
 produced under U.S. Government contract DE-AC52-06NA25396 for 
 Los Alamos National Laboratory (LANL), which is operated by 
 Los Alamos National Security, LLC for the U.S. Department of Energy. 
 The U.S. Government has rights to use, reproduce, and distribute this software.  
 NEITHER THE GOVERNMENT NOR LOS ALAMOS NATIONAL SECURITY, LLC MAKES ANY 
 WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF 
 THIS SOFTWARE.  If software is modified to produce derivative works, 
 such modified software should be clearly marked, so as not to confuse 
 it with the version available from LANL.
 
 
 Additionally, this program is free software; you can redistribute it 
 and/or modify it under the terms of the GNU General Public License 
 as published by the Free Software Foundation; version 2.0 of the License. 
 Accordingly, this program is distributed in the hope that it will be useful, 
 but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY 
 or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License 
 for more details
*/

#ifndef _CHI2DRIVER_H_
#define _CHI2DRIVER_H_

#define MAXPOINTS 100000
#define ECHOSIZE 256

extern int ggsl_errno;      /* indicates that an error has been handled */
extern FILE* errfile;

extern char* version_chi2driver;

typedef enum {eOutfile, eOutadjfile, eErrfile, eDiagnosticfile, nOutputfiles} FNUM;


int main(int argc, char *argv[]);

void    printBestChi(FILE* outfile, int fAdjusted,
                const double *sampletimes, const float *vals, const float *stderrs, 
                long ndata, int nharmonics,
                float chi0, float chid, float detrendcoeffs[], int detrendorder,
                double fmin, double deltaf, double t0, float chimargin,
                float *chireductionvalues, int nchivalues,
                char *echostring);

void    printDiagnosticFile(FILE* outfile, long ndata,
                            float chi0, float chid,
                            double freqmin, double deltaf,
                            float *chireductionvalues, int nchivalues,
                            char *echostring);


/* Read the next source in the infile, and return the number of */
/* measurements (0 for err or eof) */
int     readsource(FILE* infile, 
                char *echo, double *sampletimes, float*vals, float *stderrs);

/* Find the chisquared value for deviation from the mean */
/* if pmean == NULL, then deviation from 0 */
double chisqconst(int ndata, float *values, float *stderrs, float* pmean);


double calcCurve(float phase, int nharmonics, float* harmweights);


/* Slowly calculated harmonic coefficients at a frequency, return chisquared */
double slowharmonics(double t0, int nharmonics, float frequency, int ndata,
            double *sampletimes, float *vals, float *stderrs, float *modelvalues,
            float* coeffs);

/* Take data and subtract the best-fit polynomial, return chisquared */
double detrend(double t0, int order, int ndata,
            double *sampletimes, float *vals, float *stderrs,
            float *adjusted, float* coeffs);  

void printheader(FILE* outfile, int argc, char** argv, int detrendorder, int nharmonics, FNUM filenum);

void usage(char* name);

void my_gsl_error_handler(const char* reason, const char* file, int line, int gsl_errno);

/* choose a format specifier for frequencies that has specific precision and range */
void SetFreqFormat(double df, double freqmax); 
extern char freqformat[]; /* initialized to default frequency format */

#endif
