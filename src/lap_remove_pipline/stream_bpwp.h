
#ifndef _STREAM_BPWP_H_
#define _STREAM_BPWP_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_linalg.h>

#define MAX_POINTS 10000

typedef struct {
    double timestamp;
    double value;
    double stderr;
} DataPoint;

typedef struct {
    DataPoint buffer[MAX_POINTS];
    int head;
    int count;
} DataStream;

void init_stream(DataStream* stream);
void add_point(DataStream* stream, double timestamp, double value, double stderr);
void process_stream(DataStream* stream, FILE* output);

#endif
