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

void write_metadata_pt(const char* sim_name, int N, int M, double* a_array, double* beta_array, int n_measures) {
    char filename[512];
    snprintf(filename, sizeof(filename), "data/%s/metadata.csv", sim_name);

    FILE* file = fopen(filename, "w");
    if (file == NULL) {
        perror("Errore nella creazione di metadata.csv per il Parallel Tempering");
        return;
    }

    // Intestazione compatibile con analyze_3.py (che cerca N, a, beta)
    fprintf(file, "N,a,beta,n_measures\n");

    // Genera una riga per ogni singola replica PT
    for (int k = 0; k < M; k++) {
        fprintf(file, "%d,%.10f,%.10f,%d\n", N, a_array[k], beta_array[k], n_measures);
    }

    fclose(file);
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
