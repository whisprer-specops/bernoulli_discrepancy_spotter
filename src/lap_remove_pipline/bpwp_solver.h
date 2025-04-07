
#ifndef _BPWP_SOLVER_H_
#define _BPWP_SOLVER_H_

#include "stream_bpwp.h"
#include <gsl/gsl_vector.h>

void basis_pursuit_denoise(DataStream* stream, gsl_vector* frequencies, gsl_vector* amplitudes, int max_freqs);

#endif
