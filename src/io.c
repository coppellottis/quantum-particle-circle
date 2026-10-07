#include <stdio.h>
#include <stdlib.h>
#include "io.h"

void write_metadata(const char* sim_name, int phase, double fixed_val, 
                    const int* N_array, int n_N, int n_measures) {
    
    char metadata_path[512];
    snprintf(metadata_path, sizeof(metadata_path), "data/%s/metadata.csv", sim_name);

    FILE* metadata = fopen(metadata_path, "w");
    if (metadata == NULL) {
        perror("WARNING: Non è stato possibile aprire il file delle misure in scrittura.");
        exit(EXIT_FAILURE);
    }

    // Header del CSV 
    fprintf(metadata, "phase,N,a,beta,n_measures\n");
    
    for (int i = 0; i < n_N; i++) {
        int N = N_array[i];
        double a, beta;

        // Calcoliamo i parametri dipendenti in base alla fase scelta
        if (phase == 1) {
            a = fixed_val;
            beta = N * a;
        } else {
            beta = fixed_val;
            a = beta / (double)N;
        }

        // Scrittura della riga
        fprintf(metadata, "%d,%d,%f,%f,%d\n", phase, N, a, beta, n_measures);
    }

    fclose(metadata);

    return;
}

FILE* open_measurement_file(const char* sim_name, int N, double a, double beta) {
    char filename[256];
    snprintf(filename, sizeof(filename), "data/%s/N%d_a%.4f_beta%.4f.csv", sim_name, N, a, beta);

    FILE* file = fopen(filename, "w");
    if (file == NULL) {
        perror("WARNING: Non è stato possibile aprire il file delle misure in scrittura.");
        exit(EXIT_FAILURE);
    }

    fprintf(file, "sweep,E,Q\n");
    return file;
}
