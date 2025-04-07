#include "fastchi2.h"  // Assuming fastchi2.h defines chisqconst and other required functions
#include <stdio.h>
#include <stdlib.h>

#define WINDOW_SIZE 100  // Sliding window size for analysis
#define ALARM_THRESHOLD 3.84  // Chi-squared threshold for 1 degree of freedom (95% confidence level)

// Pre-allocated buffers
float *results_window;
float *stderrs;
float *adjvalues;
int total_trials = 0;
int total_wins = 0;

// Initialize buffers for sliding window analysis
void initialize_buffers() {
    results_window = (float *)malloc(WINDOW_SIZE * sizeof(float));
    stderrs = (float *)malloc(WINDOW_SIZE * sizeof(float));
    adjvalues = (float *)malloc(WINDOW_SIZE * sizeof(float));

    if (!results_window || !stderrs || !adjvalues) {
        fprintf(stderr, "Memory allocation error.\n");
        exit(1);
    }

    // Initialize standard errors (constant for Bernoulli trials)
    for (int i = 0; i < WINDOW_SIZE; i++) {
        stderrs[i] = 1.0f;
    }
}

// Process each incoming packet outcome (1 for win, 0 for loss)
void process_packet(int outcome) {
    static int current_index = 0;

    // Update sliding window with the new outcome
    results_window[current_index] = outcome == 1 ? 1.0f : 0.0f;  // Win = 1, Loss = 0
    current_index = (current_index + 1) % WINDOW_SIZE;  // Circular buffer index
    total_trials++;

    // Perform chi-squared analysis when enough data is available
    if (total_trials >= WINDOW_SIZE) {
        float mean0;
        float chi0 = chisqconst(WINDOW_SIZE, results_window, stderrs, &mean0);

        if (chi0 > ALARM_THRESHOLD) {
            printf("ALARM: Chi-squared deviation detected! Chi: %f\n", chi0);
        } else {
            printf("Stable: Chi: %f\n", chi0);
        }
    }
}

// Simulate packet stream processing (replace this with real HTTP packet handling)
void simulate_http_stream() {
    // Dummy data stream (1 for win, 0 for loss)
    int outcome_stream[] = {1, 0, 1, 1, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 1, 1, 0, 1, 1, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 0, 1, 1, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 1};, 
    int stream_size = sizeof(outcome_stream) / sizeof(outcome_stream[0]);

    // Process each packet (outcome)
    for (int i = 0; i < stream_size; i++) {
        process_packet(outcome_stream[i]);
    }
}

int main() {
    // Initialize buffers for sliding window and chi-squared analysis
    initialize_buffers();

    // Simulate HTTP packet stream processing
    simulate_http_stream();

    // Free allocated memory buffers
    free(results_window);
    free(stderrs);
    free(adjvalues);

    return 0;
}