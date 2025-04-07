
#include "stream_bpwp.h"

void init_stream(DataStream* stream) {
    stream->head = 0;
    stream->count = 0;
}

void add_point(DataStream* stream, double timestamp, double value, double stderr) {
    int idx = (stream->head + stream->count) % MAX_POINTS;
    stream->buffer[idx].timestamp = timestamp;
    stream->buffer[idx].value = value;
    stream->buffer[idx].stderr = stderr;
    if (stream->count < MAX_POINTS) stream->count++;
    else stream->head = (stream->head + 1) % MAX_POINTS;
}

void process_stream(DataStream* stream, FILE* output) {
    gsl_matrix* X = gsl_matrix_alloc(stream->count, 2);
    gsl_vector* y = gsl_vector_alloc(stream->count);
    gsl_vector* c = gsl_vector_alloc(2);
    gsl_matrix* cov = gsl_matrix_alloc(2, 2);

    for (int i = 0; i < stream->count; i++) {
        double t = stream->buffer[(stream->head + i) % MAX_POINTS].timestamp;
        gsl_matrix_set(X, i, 0, 1.0);
        gsl_matrix_set(X, i, 1, t);
        gsl_vector_set(y, i, stream->buffer[(stream->head + i) % MAX_POINTS].value);
    }
    
    gsl_multifit_linear_workspace* work = gsl_multifit_linear_alloc(stream->count, 2);
    gsl_multifit_linear(X, y, c, cov, NULL, work);

    double c0 = gsl_vector_get(c, 0);
    double c1 = gsl_vector_get(c, 1);
    fprintf(output, "Trend intercept: %f, slope: %f\n", c0, c1);

    // For outputting laplace noise mask instead of masked signal:
//    fprintf(output, "Timestamp\tResidual (Laplace Noise)\n");
//    for (int i = 0; i < stream->count; i++) {
//        double t = stream->buffer[(stream->head + i) % MAX_POINTS].timestamp;
//        double original_val = stream->buffer[(stream->head + i) % MAX_POINTS].value;
//        double predicted_val = c0 + c1 * t;
//        double residual = original_val - predicted_val;
//        fprintf(output, "%.6f\t%.6f\n", t, residual);
//    }

    gsl_multifit_linear_free(work);
    gsl_matrix_free(X);
    gsl_vector_free(y);
    gsl_vector_free(c);
    gsl_matrix_free(cov);
}
