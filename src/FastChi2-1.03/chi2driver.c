/*
* fastchi driver program
* main() for running fast-chi algorithm described in
* Palmer, D.M., 2007 (submitted)
* "A Fast Chi-squared Technique For Period Search of Irregularly Sampled Data."
*
* Usage information follows the copyright notice.
*/
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
/* Arguments
    runchi2 nharmonics freqmax
            [-f freqmin]
            [-d detrendorder]
            [-t t0]
            [-T timespan]
            [-s oversampling]
            [-m chimargin]
            [-i infile]
            [-o outfile]
            [-O outfile_with_frequency_adjustment]
            [-e stderrfile]
            [-D diagnosticfile]
            [-c 'comment string']
    nharmonics includes fundamental. ( 3 -> {fundamental, 2f, 3f})
    freqmax is the highest fundamental frequency searched and determines FFT size
    freqmin limits the search range (defaults to 1 cycle over the span of each observation)
    detrendorder is the order of the polynomial detrending, default 3,
        value of zero will just do constant
    t0 is the time zero of the sinusoids and the detrending polynomial, 
        default is midway between endpoints.  Value is given in each output line
    timespan is duration used to determine order of the FFT.  If it is not given, it is taken
        from the span of the first source
    oversampling determines how tight the search is. The fundamental frequency step
         is 1/(oversampling * nharmonics * firsttimespan), rounded down
         to give a power of 2 size in the FFTs
    chimargin is how much improvement to expect from the singleChi peaking 
         routine.
    infile format is discussed in next comment block
    outfile and outfile_with_frequency_adjust may both be requested
        outfile_with_frequency_adjust finds the frequency with the greatest significance
        after adjusting for number of trials to that point (taken as proportional to the frequency),
        as a Frequentist would find appropriate.
        If you are a Bayesian, it uses the likelihood based on a prior probability function
        that is uniform in log(frequency) (== uniform in log(period)), and this works out to exactly
        the same calculation.
        If neither is requested, then outfile is set to stdout.  Stdout can be explicitly requested
        with name '-'.
        Output file format is discussed two comment blocks down
    diagnosticfile contains the chi-squared reduction (larger is better) for each frequency
    The comment string is ignored.  However since the entire command line is included
        in the outfile header, it is echoed there.
*/
/* 
The input file can include any number of sources, with a set of lines for each source in the format:
    echoline
    n
    t0 x0 err0
    t1 x1 err1
    ...
    tn xn errn
The echoline is used as the first field of the corresponding line of the output file, and so for ease of
interpretation it should not have any tab characters in it and it should not start with '#'.
*/

/* 
The output file format is a multi-line header (lines beginning with #)
followed by a tab-delimited line for each source in the format:
    echoline fbest ndata chimin chi0 chidetrend tzero harmbaseline sin1 cos1 sin2 cos2...  d0 d1 d2...
where 
echoline comes from the input file (e.g. the name of the source)
fbest is the best frequency, 
ndata is the number of data samples that went into the calculation,
chimin is the minimum chisquared at fbest, 
chidetrend is the chisquared of the detrended data,
chi0 is the chisquared of the non-detrended data
    (== chidetrend if detrendorder == 0)
tzero is the (arbitrary) phase=0 time
harmbaseline is the constant term
sin1 cos1 are the Fourier components at fbest 
sin2 cos2 are the Fourier components at 2*fbest
...
d0 d1 d2... coefficients of detrending polynomial fit
Note that there are two constant terms (harmbaseline and d0) that must be added.

As an example, if you wish to calculate the model value at time t, using
parameters named as in the header line, for a 3-harmonic fit 
with an order-2 polynomial detrending:

phi(t) = (t - tzero) * fbest ; // in cycles

Fit(t) = (harmbaseline + d0)
        + sin1 * sin(1 * phi(t) * 2 * pi) + cos1 * cos(1 * phi(t) * 2 * pi)
        + sin2 * sin(2 * phi(t) * 2 * pi) + cos2 * cos(2 * phi(t) * 2 * pi)
        + sin3 * sin(3 * phi(t) * 2 * pi) + cos3 * cos(3 * phi(t) * 2 * pi)
        + d1 * (t - tzero) + d2 * (t - tzero) * (t - tzero);
*/


#include <unistd.h>
#include <stdio.h>
#include "fastchi2.h"
#include "chi2driver.h"
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_errno.h>
#include <assert.h>

char* version_chi2driver = "1.03";

int ggsl_errno = 0;
FILE* errfile;


int main(int argc, char *argv[])
{
/* These are allocated static because some compilers do not like large arrays */
/* allocated on the stack.  FIXME: They should be malloc()ed with MAXPOINTS   */
/* allowed as a command line argument, or they should be grown as required    */
/* based on the number of points in the actual input.  The times are double   */
/* and the values are floats because that is typically the appropriate        */
/* precision of these kinds of measurements, especially when the times are    */
/* relative to a long-past epoch */
    static double  sampletimes[MAXPOINTS];
    static float   values[MAXPOINTS];
    static float   stderrs[MAXPOINTS];
    static float   adjvalues[MAXPOINTS];   /* values after detrending */
    long idatum;
  
    char echostring[ECHOSIZE+1];

    int nharmonics;
    double freqmax;
    double freqmin = 0;
    int detrendorder = 0;
    double t0 = -1;
    float oversampling = 1.0;
    float obsrange = 0.0;   /* duration of observation */
    double tstart;           /* time at which to start FFT arrays (half a bin before earliest datum) */
    double tmindata;        /* time of first datum */
    double tmaxdata;        /* time of last datum */
  
    float chimargin = 0.0;  /* How much improvement in Chi^2 to allow for by singleChi() */
    FILE* infile = stdin;
    FILE* outfile = NULL;
    FILE* outadjfile = NULL;
    FILE* diagnosticfile = NULL;
    FILE*  allOutFiles[nOutputfiles];
    FNUM fnum;

    long    nrealpoints = 0;
    long    nchivalues = 0;


    float   *data_var = NULL;  /* [nrealpoints] of data/variance */
    float   *one_var = NULL;   /* [nrealpoints] of 1/variance */
    float   *chireductionvalues = NULL; /* [nrealpoints/(2 * nharmonics)] */

    double deltaf;
    double deltat;
    
    int ndata;

    int t0given = 0;    /* set for explicit t0 value */


    gsl_set_error_handler(my_gsl_error_handler);

    errfile = stderr;
    
    if (argc < 3)
    {
        usage(argv[0]);
    }
    else
    {
        char ch;
        if (1 != sscanf(argv[1], "%d", &nharmonics) || 1 != sscanf(argv[2], "%lf", &freqmax))
        {
            usage(argv[0]);
            exit(-1);
        }
        while ((ch = getopt(argc - 2, argv + 2, "f:d:t:T:s:m:pi:o:O:e:D:c:")) != -1)
        {
            switch (ch)
            {
                case 'f': 
                    if (1 != sscanf(optarg, "%lf", &freqmin))
                        usage(argv[0]);
                    break;
                case 'd': 
                    if (1 != sscanf(optarg, "%d", &detrendorder))
                        usage(argv[0]);
                    break;
                case 't': 
                    if (1 != sscanf(optarg, "%lf", &t0))
                        usage(argv[0]);
                    t0given = 1;
                    break;
                case 'T':
                    if (1 != sscanf(optarg, "%f", &obsrange))
                        usage(argv[0]);
                    break;
                case 's':
                    if (1 != sscanf(optarg, "%f", &oversampling))
                        usage(argv[0]);
                    break;
                case 'i': 
                    if (0 == (infile = fopen(optarg, "rb")))
                        usage(argv[0]);
                    break;
                case 'm':
                    if (1 != sscanf(optarg, "%f", &chimargin))
                        usage(argv[0]);
                    break;
                case 'o':
                    if (0 == strcmp(optarg, "-"))
                    {
                        outfile = stdout;
                    }
                    else if (0 == (outfile = fopen(optarg, "wb")))
                        usage(argv[0]);
                    break;
                case 'O': 
                    if (0 == strcmp(optarg, "-"))
                    {
                        outadjfile = stdout;
                    }
                    else if (0 == (outadjfile = fopen(optarg, "wb")))
                        usage(argv[0]);
                    break;
                case 'e': 
                    if (0 == (errfile = fopen(optarg, "wb")))
                        usage(argv[0]);
                    break;
                case 'D':
                    if (0 == (diagnosticfile = fopen(optarg, "wb")))
                        usage(argv[0]);
                    break;
                case 'c':
                    break;
                case '?':
                default:
                    usage(argv[0]);
            }
        }
    }
    
    if (outfile == NULL && outadjfile == NULL)
    {   /* if no output specified, default to send non-frequency-adjusted output to stdout */
        outfile = stdout;     
    }
    
    allOutFiles[eOutfile] = outfile;
    allOutFiles[eOutadjfile] = outadjfile;
    allOutFiles[eErrfile] = errfile;
    allOutFiles[eDiagnosticfile] = diagnosticfile;
    
    
    for (fnum = eOutfile ; fnum < nOutputfiles ; fnum++)
    {
        if (allOutFiles[fnum])
        {
            printheader(allOutFiles[fnum], argc, argv, detrendorder, nharmonics, fnum);
        }        
    }
    

    /* Read the first observation with 2 or more data points in it */
    while (-1 != (ndata = readsource(infile, echostring, sampletimes, values, stderrs))
          && 2 > ndata)
    {
        fprintf(errfile, "%s: Insufficient data (%d points)\n", echostring, ndata);
    }
    
    if (ndata < 2) usage(argv[0]);  /* quit if none of the data in the file have 2 or more points */

  
    /* The obsrange, span of observing times can be passed as a -T argument, 
    /* or it can be taken from the first source processed in a run. */
    if (obsrange == 0)
    {
        /* Find the earliest and latesttimes of the data */
        tmindata = sampletimes[0];
        tmaxdata = sampletimes[0];
        for (idatum = 0 ; idatum < ndata ; idatum++)
        {
          if (tmindata > sampletimes[idatum]) tmindata = sampletimes[idatum]; 
          if (tmaxdata < sampletimes[idatum]) tmaxdata = sampletimes[idatum]; 
        }
        obsrange = tmaxdata - tmindata;
    }


    {   /* IMPROVEME: put this in a separate function */
            /* returns deltat, deltaf, data_var, one_var, chireductionvalues */
        float mintotalrange = obsrange * nharmonics * oversampling * 2;  
                /* the timespan of the data and the padding to higher frequencies */
        float totalrange;
        
        deltat = 1/(2*freqmax * 2*nharmonics);
            /* Nyquist spacing for freqmax * 2*nharamonics .  2*nharmonics is required to get */
            /* cross terms like cos(nf)*cos(nf) which includes cos(2nf) etc. */
        
        nrealpoints = 1L << (int)(ceil(log2(mintotalrange/deltat)));
        if (nrealpoints <= 1)
        {
          fprintf(errfile, "Too many time bins (%.0f) in FFT calculation\n", mintotalrange/deltat);
          assert(0);  /* Fails if too many time bins to be indexed by a long int: */
                          /* when mintotalrange/deltat >= 2^30 ~ 1e9 if long int is int_32 */
        }
      
        totalrange = nrealpoints * deltat;
        deltaf = 1/totalrange;
        
        SetFreqFormat(deltaf, freqmax);
      
        for (fnum = eOutfile ; fnum < nOutputfiles ; fnum++) 
        {
            if (allOutFiles[fnum])
            {
                fprintf(allOutFiles[fnum], "# %ld real timeslots of %g time-units,\n# (%ld / %.g time-units of observation)\n",
                        nrealpoints, deltat, (long)(nrealpoints * obsrange/totalrange), obsrange);
                fprintf(allOutFiles[fnum], "# giving a %g cycles/time-unit (%g/observation) freq. grid\n",
                        deltaf, deltaf * obsrange);
            }
        }
        
        nchivalues = nrealpoints / (2 * 2 * nharmonics);
        
        data_var = malloc(nrealpoints * sizeof(float));
        one_var = malloc(nrealpoints * sizeof(float));
        chireductionvalues = malloc(nchivalues * sizeof(float));
        if (!data_var || !one_var || !chireductionvalues)
        {
            fprintf(errfile, "Unable to allocate sufficient memory\n");
            assert(0);
        }
    }

    do  /* Over all datasets */
    {
        if (ndata > (detrendorder + 2 * nharmonics + 2))
        {
            float chi0; /* chisquared of constant value */
            float chid; /* chisquared of detrended data */
            float mean0;
            
            float   detrendcoeffs[MAXORDER+1];
            
            /* Find the earliest and latesttimes of the data */
            tmindata = sampletimes[0];
            tmaxdata = sampletimes[0];
            for (idatum = 0 ; idatum < ndata ; idatum++)
            {
              if (tmindata > sampletimes[idatum]) tmindata = sampletimes[idatum]; 
              if (tmaxdata < sampletimes[idatum]) tmaxdata = sampletimes[idatum]; 
            }
          
            if (!t0given)
            {
                /* Choose a t0 that is the midway between the earliest and latest data */
                t0 = (tmindata + tmaxdata)/2;
            }
               
            tstart = tmindata - deltat/2; /* First datum is in the middle of bin[0] */

            chi0 = chisqconst(ndata, values, stderrs, &mean0);

            chid = detrend(t0, detrendorder, ndata, sampletimes, values, stderrs, adjvalues, detrendcoeffs);                    
            
            
            fastChi(sampletimes, adjvalues, stderrs,
                    ndata, nharmonics, deltat,
                    nrealpoints, nchivalues, tstart,
                    data_var, one_var, chireductionvalues);
            
            /* for (fnum = 0 ; fnum <= 1 ; fnum++) */    /* outfile and adjoutfile */
            for (fnum = eOutadjfile ; fnum >= eOutfile && fnum < nOutputfiles ; fnum--)
            {                       /* This double test is required because enums may be unsigned */
                if (allOutFiles[fnum])
                {
                    printBestChi(allOutFiles[fnum], (eOutadjfile == fnum),
                            sampletimes, adjvalues, stderrs,
                            ndata, nharmonics,
                            chi0, chid, detrendcoeffs, detrendorder,
                            freqmin, deltaf, t0, chimargin,
                            chireductionvalues, nchivalues,
                            echostring);
                }
            }
            
            if (allOutFiles[eDiagnosticfile])
            {
                printDiagnosticFile(allOutFiles[eDiagnosticfile],
                                    ndata, chi0, chid,
                                    freqmin, deltaf,
                                    chireductionvalues, nchivalues,
                                    echostring);                
            }
        }
        else
        {
            fprintf(errfile, "%s: Insufficient data (%d points)\n", echostring, ndata);
            fflush(errfile);
        }
    } while (-1 != (ndata = readsource(infile, echostring, sampletimes, values, stderrs)));

    fprintf(errfile, "# Done\n");
    for (fnum = 0 ; fnum < nOutputfiles ; fnum++) 
    {
        if (allOutFiles[fnum])
        {
            fclose(allOutFiles[fnum]);
        }        
    }
    return 0;
}


/* Find and print the parameters with the best chisquared, or, if fAdjusted is true, */
/* the best trial-adjusted probability: */
/* smallest P(> \DELTA\chi^2, \nu = 2*nharmonics) * freq */

void    printBestChi(FILE* outfile, int fAdjusted,
                const double *sampletimes, const float *vals, const float *stderrs, 
                long ndata, int nharmonics,
                float chi0, float chid, float detrendcoeffs[], int detrendorder,
                double freqmin, double deltaf, double t0, float chimargin,
                float *chireductionvalues, int nchivalues,
                char *echostring)
{
    float   harmweights[2*MAXORDER + 1];
    float  best_peakup;      /* Maximum improvement of chisquared from searching off the gridpoints */ 
    static double best_peakup_across_stars = 0; /* Best peakup for any star this run */
    
    double freqbest;
    
    freqbest = findBestFrequency(fAdjusted,
                                 sampletimes, vals, stderrs, 
                                 ndata, nharmonics,
                                 freqmin,  deltaf,  t0,  chimargin,
                                 chireductionvalues,  nchivalues, 
                                 &best_peakup);

    if (best_peakup > best_peakup_across_stars)
    {
        fprintf(errfile, "# %s : Ungridded search peaked-up chi^2 by %0.2f . %s\n",
                echostring,
                best_peakup,
                (best_peakup > chimargin ? "*** THIS EXCEEDS -m CHI-MARGIN ARGUMENT ***" : ""));
        fflush(errfile);
        best_peakup_across_stars = best_peakup;
    }
    
    if (ggsl_errno)
    {
        ggsl_errno = 0;
        fprintf(errfile, "%s: No minimum found due to error.\n", echostring);
        fflush(errfile);
    }
    else
    {
        int i;
        double chireduction = singleChi(sampletimes, vals, stderrs, 
                                ndata, nharmonics, t0,
                                freqbest, harmweights);
        fprintf(outfile, "%s\t", echostring);
        fprintf(outfile, freqformat, freqbest);
        fprintf(outfile, "\t%d\t%f\t%f\t%f\t%f",
                            (int)ndata, chid - chireduction, chi0, chid, t0);
        for (i = 0 ; i <= 2*nharmonics ; i++)
            fprintf(outfile, "\t%g", harmweights[i]);
        for (i = 0 ; i <=  detrendorder ; i++)
            fprintf(outfile, "\t%g", detrendcoeffs[i]);
        fprintf(outfile, "\n");
    }   
    fflush(outfile);
}

/* Find and print the parameters with the best chisquared, or, if fAdjusted is true, */
/* the best trial-adjusted probability: */
/* smallest P(> \DELTA\chi^2, \nu = 2*nharmonics) * freq */

void    printDiagnosticFile(FILE* outfile, long ndata,
                     float chi0, float chid,
                     double freqmin, double deltaf,
                     float *chireductionvalues, int nchivalues,
                     char *echostring)
{
    int fi;
    int fimin = freqmin/deltaf;
    fprintf(outfile, "%s\n", echostring);
    fprintf(outfile, "%d %g %g %d\n", nchivalues - fimin, chi0, chid, (int)ndata);
    for (fi = fimin ; fi < nchivalues ; fi++)
    {
        fprintf(outfile, freqformat,  fi*deltaf);
        fprintf(outfile, "\t%g\n", chireductionvalues[fi]);
    }
    fflush(outfile);
}


/* Read the next source in the infile, and return the number */
/* measurements (0 for err or eof) */
int     readsource(FILE* infile, 
                char *echo, double *sampletimes, float *vals, float *stderrs)
{
    int i;
    int nmeas = 0;
        /* 256 == ECHOSIZE */
    fscanf(infile, "%256[\n\r]", echo);
    if (1 == fscanf(infile, "%256[^\n\r]", echo)
        && 1 == fscanf(infile, "%d", &nmeas)
        && nmeas <= MAXPOINTS)
    {
        for (   i = 0 
              ; i < nmeas 
                && 3 == fscanf(infile, "%lf %f %f", 
                                        &sampletimes[i], &vals[i], &stderrs[i])
              ; i++)
        {
            ;  /* read done in for loop header */
        }
        if (i == nmeas)
            return nmeas;
    }
  
    assert(nmeas <= MAXPOINTS); /* Throw an assert if that was the problem */
                                /* IMPROVEME: Make number of points dynamic */
    return -1;   /* either end of file or something wrong */
}


           
double chisqconst(int ndata, float *values, float *stderrs, float* pmean)
{                    /* not worth optimizing away the obvious inefficiencies */
    float mean;
    float chisq = 0;
    int i;
    if (pmean)
    {
        float sumvalovar = 0;
        float sum1ovar = 0;
        for (i = 0 ; i < ndata ; i++)
        {
            float invvar = 1/(stderrs[i]*stderrs[i]);
            sum1ovar += invvar;
            sumvalovar += values[i]*invvar;
        }
        mean = sumvalovar / sum1ovar;
        *pmean = mean;
    }
    else
    {
        mean = 0;
    }
    
    for (i = 0 ; i < ndata ; i++)
        chisq += (values[i]-mean)*(values[i]-mean)/(stderrs[i]*stderrs[i]);
        
    return chisq;
}


/* Slowly calculated harmonic coefficients at a frequency, return chisquared */
double slowharmonics(double t0, int nharmonics, float frequency , int ndata,
            double *sampletimes, float *vals, float *stderrs, float *modelvalues,
            float* coeffs)
{
    int i, j;
    float omega = frequency * 2 * M_PI;
    float chisqactual = 0;
    gsl_multifit_linear_workspace   *pwork = gsl_multifit_linear_alloc(ndata, 2*nharmonics+1);
    gsl_matrix  *x = gsl_matrix_alloc(ndata, 2*nharmonics+1);
    gsl_vector  *y = gsl_vector_alloc(ndata);
    gsl_vector  *c = gsl_vector_alloc(2*nharmonics+1);
    gsl_matrix  *cov = gsl_matrix_alloc(2*nharmonics+1, 2*nharmonics+1);
    

    double  chisq;

    for (i = 0 ; i < ndata ; i++)
    {
        double dt = sampletimes[i] - t0; 
        double phi = dt * omega; 
        for (j = 0 ; j <= 2*nharmonics ; j++)
        {
            gsl_matrix_set(x, i, j, Mfunc(j, phi)/stderrs[i]);
        }
        gsl_vector_set(y, i, vals[i]/stderrs[i]);
    }

    gsl_multifit_linear(x, y, c, cov, &chisq, pwork);
    
            /* move the coefficients to someplace usable */
    for (j = 0 ; j <= 2*nharmonics ; j++)
    {
        coeffs[j] = gsl_vector_get(c, j);
    }
    
    for (i = 0 ; i < ndata ; i++)
    {
        modelvalues[i] = 0;
        for (j = 0 ; j <= 2*nharmonics ; j++)
        {
            modelvalues[i] += coeffs[j] * gsl_matrix_get(x,i,j);
        }
        modelvalues[i] *= stderrs[i];
        chisqactual += (modelvalues[i]-vals[i])*(modelvalues[i]-vals[i])/(stderrs[i]*stderrs[i]);
    }
            
    gsl_multifit_linear_free(pwork);
    gsl_matrix_free(x);
    gsl_vector_free(y);
    gsl_vector_free(c);
    gsl_matrix_free(cov);
    
    return chisq;
}

/* Take data and subtract the best-fit polynomial, return chisquared */
double detrend(double t0, int order, int ndata,
            double *sampletimes, float *vals, float *stderrs,
            float *adjusted, float* coeffs)
{
    int i, j;
    double sumdov = 0;
    double sum1ov = 0;
    double  chisq = 0;
    float ave;
    
    for (i = 0 ; i < ndata ; i++)
    {
        double ivar = 1/stderrs[i];
        ivar *= ivar;
        sumdov += vals[i] * ivar;
        sum1ov += ivar;
    }
    ave = sumdov/sum1ov;
    
    if (order == 0)     /* average only */
    {
        coeffs[0] = ave;
        for (i = 0 ; i < ndata ; i++)
        {
            adjusted[i] = vals[i] - ave;
            chisq += (adjusted[i] * adjusted[i])/(stderrs[i] * stderrs[i]);
        }
        return chisq;
    }
    else
    {
        gsl_multifit_linear_workspace   *pwork = gsl_multifit_linear_alloc(ndata, order+1);
        gsl_matrix  *x = gsl_matrix_alloc(ndata, order+1);
        gsl_vector  *y = gsl_vector_alloc(ndata);
        gsl_vector  *c = gsl_vector_alloc(order+1);
        gsl_matrix  *cov = gsl_matrix_alloc(order+1, order+1);


        double  chisqactual;

        for (i = 0 ; i < ndata ; i++)
        {
            double v = 1/stderrs[i];    /* rather than use weight matrix, measure in sigmas */
            double dt = sampletimes[i] - t0; 
            for (j = 0 ; j <= order ; j++)
            {
                gsl_matrix_set(x, i, j, v);
                v *= dt; 
            }
            gsl_vector_set(y, i, vals[i]/stderrs[i]);
        }

        gsl_multifit_linear(x, y, c, cov, &chisq, pwork);
        
                /* move the coefficients to someplace usable */
        for (j = 0 ; j <= order ; j++)
        {
            coeffs[j] = gsl_vector_get(c, j);
        }
        
        chisqactual = 0;
        for (i = 0 ; i < ndata ; i++)
        {
            float v = coeffs[order];
            float dt = sampletimes[i] - t0;
            for (j = order-1 ; j >= 0 ; j--)
            {
                v = coeffs[j] + v * dt; /* evaluate polynomial */
            }
            
            adjusted[i] = vals[i] - v;  /* adjust by subtracting polynomial estimate */
            chisqactual += (adjusted[i] * adjusted[i])/(stderrs[i]*stderrs[i]);
        }
        
        gsl_multifit_linear_free(pwork);
        gsl_matrix_free(x);
        gsl_vector_free(y);
        gsl_vector_free(c);
        gsl_matrix_free(cov);
        
        return chisq;
    }
}

double calcCurve(float phase, int nharmonics, float* harmweights)
{
    int i;
    double v = harmweights[0];
    for (i = 1 ; i <= nharmonics ; i++)
    {
        v += harmweights[2*i - 1]*sin(2*M_PI*phase*i)
            + harmweights[2*i]*cos(2*M_PI*phase*i);
    }
    return v;
}

void printheader(FILE* outfile, int argc, char** argv, int detrendorder, int nharmonics, FNUM filenum)
{
    int i;
    fprintf(outfile, "#");
    for (i = 0 ; i < argc ; i++)
    {
        fprintf(outfile, strchr(argv[i], ' ') ? " '%s'" : " %s", argv[i]);
    }
    switch (filenum) 
    {
        case eOutfile :
            fprintf(outfile, "\n# Unadjusted chi-squared minimum");
            break;
        case eOutadjfile :
            fprintf(outfile, "\n# Frequency-adjusted chi-squared minimum");
            break;
        case eErrfile :
            fprintf(outfile, "\n# Error output");
            break;
        case eDiagnosticfile :
            fprintf(outfile, "\n# Diagnostic output");
            break;
        case nOutputfiles :
        default :
            assert(0);
            break;
    }
    if (filenum <= eOutadjfile)
    {
        fprintf(outfile, "\n# Columns:echostring\tfbest\tndata\tchimin\tchi0\tchidetrend\ttzero\tharmbaseline");        
        
        for (i = 1 ; i <= nharmonics ; i++)
        {
            fprintf(outfile, "\tsin%d\tcos%d", i, i);
        }
        for (i = 0 ; i <=  detrendorder ; i++)
        {
            fprintf(outfile, "\td%d", i);
        }
    }    
    fprintf(outfile, "\n# %s\n# (fastchi2.c version %s, chi2driver.c version %s, compiled %s)\n", 
                            argv[0], version_fastchi2, version_chi2driver, __DATE__);
    fflush(outfile);
}

void usage(char* name)
{
    fprintf(errfile, "Usage: %s nharmonics freqmax\n"
                    "\t[-f freqmin]\n"
                    "\t[-d detrendorder]\n"
                    "\t[-t t0]\n"
                    "\t[-T timespan]\n"
                    "\t[-s oversampling]\n"
                    "\t[-m chi2margin]\n"
                    "\t[-i infile]\n"
                    "\t[-o outfile]\n"
                    "\t[-O outfile_with_frequency_adjustment]\n"
                    "\t[-e stderrfile]\n"
                    "\t[-D diagnosticfile]\n"
                    "\t[-c 'comment string']\n",
             name);
    exit(-1);
}


void my_gsl_error_handler(const char* reason, const char* file, int line, int gsl_errno)
{
    if (gsl_errno != GSL_EMAXITER)  /* ignore this if due to Chi-squared cdf out of range */
    {
        fprintf(errfile, "GSL ERROR: %s (%d), in %s, line %d\n", reason, gsl_errno, file, line);
        ggsl_errno = gsl_errno;
    }
}


char freqformat[10] = "%g";  /* initialized to default frequency format */

void SetFreqFormat(double df, double freqmax)
{
    if ((freqmax > 1e8 && df > 1e3) || (freqmax < 1e-2))
    {
        /* Format should be X.XXXXeSxx */
        int nsigfigs = 2 + floor(log10(freqmax)) - floor(log10(df));
        snprintf(freqformat, 9, "%%.%1de", nsigfigs-1);
          /* %e precision is number of digits after decimal */
    }
    else
    {
        /* Digits before and after decimal */
        int nafter = 1 - floor(log10(df));
        int nbefore = 1 + floor(log10(freqmax));
        if (nbefore < 1) 
        {
            nbefore = 1; /* 0.xxxx */
        }
        if (nafter < 1)
        {
            snprintf(freqformat, 9, "%%%d.0f", nbefore);
        }
        else
        {
            snprintf(freqformat, 9, "%%%d.%df", nbefore+nafter+1, nafter);
        }
    }
}
