
#include "bpwp_solver.h"

void basis_pursuit_denoise(DataStream* stream, gsl_vector* frequencies, gsl_vector* amplitudes, int max_freqs) {
    size_t N = stream->count;
    gsl_matrix* A = gsl_matrix_alloc(N, max_freqs);
    gsl_vector* b = gsl_vector_alloc(N);
    gsl_vector* x = gsl_vector_alloc(max_freqs);

    for (size_t i = 0; i < N; i++) {
        double t = stream->buffer[(stream->head + i) % MAX_POINTS].timestamp;
        gsl_vector_set(b, i, stream->buffer[(stream->head + i) % MAX_POINTS].value);
        for (int j = 0; j < max_freqs; j++) {
            double freq = gsl_vector_get(frequencies, j);
            gsl_matrix_set(A, i, j, sin(2 * M_PI * freq * t));
        }
    }

    gsl_multifit_linear_workspace* workspace = gsl_multifit_linear_alloc(N, max_freqs);
    gsl_matrix* cov = gsl_matrix_alloc(max_freqs, max_freqs);
    gsl_multifit_linear(A, b, x, cov, NULL, workspace);

    for (int i = 0; i < max_freqs; i++) {
        gsl_vector_set(amplitudes, i, gsl_vector_get(x, i));
    }

    gsl_multifit_linear_free(workspace);
    gsl_matrix_free(A);
    gsl_vector_free(b);
    gsl_vector_free(x);
    gsl_matrix_free(cov);
}
