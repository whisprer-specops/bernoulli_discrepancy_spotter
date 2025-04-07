# I'll prepare a full C implementation (source files) for estimating Laplace noise parameters (location and scale).
# This complements your existing Laplace regression pipeline, enabling precise noise management.

import zipfile

# Create files for Laplace noise parameter estimation
estimator_h = """
#ifndef _LAPLACE_ESTIMATOR_H_
#define _LAPLACE_ESTIMATOR_H_

#include <stddef.h>

void estimate_laplace_params(const double* data, size_t n, double* location, double* scale);

#endif
"""

estimator_c = """
#include "laplace_estimator.h"
#include <stdlib.h>
#include <math.h>

int compare(const void* a, const void* b) {
    double diff = (*(double*)a - *(double*)b);
    return (diff > 0) - (diff < 0);
}

double median(double* data, size_t n) {
    qsort(data, n, sizeof(double), compare);
    if (n % 2 == 0)
        return (data[n / 2 - 1] + data[n / 2]) / 2;
    else
        return data[n / 2];
}

void estimate_laplace_params(const double* data, size_t n, double* location, double* scale) {
    double* temp = malloc(n * sizeof(double));
    for (size_t i = 0; i < n; i++) temp[i] = data[i];

    double loc = median(temp, n);
    *location = loc;

    for (size_t i = 0; i < n; i++) temp[i] = fabs(data[i] - loc);
    *scale = median(temp, n) / log(2.0);  // Scale estimation using MAD (Median Absolute Deviation)

    free(temp);
}
"""

# Save files to disk
base_dir_estimator = '/mnt/data/laplace_estimator'
os.makedirs(base_dir_estimator, exist_ok=True)

with open(f"{base_dir_estimator}/laplace_estimator.h", "w") as f:
    f.write(estimator_h)

with open(f"{base_dir_estimator}/laplace_estimator.c", "w") as f:
    f.write(estimator_c)

# Package files into a zip archive
zip_path_estimator = '/mnt/data/laplace_estimator_code.zip'
with zipfile.ZipFile(zip_path_estimator, 'w') as zipf:
    zipf.write(f"{base_dir_estimator}/laplace_estimator.h", "laplace_estimator.h")
    zipf.write(f"{base_dir_estimator}/laplace_estimator.c", "laplace_estimator.c")

zip_path_estimator
