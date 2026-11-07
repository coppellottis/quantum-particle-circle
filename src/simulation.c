#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "simulation.h"
#include "polymer.h"
#include "metropolis.h"
#include "io.h"
#include "progress.h"
#include "observables.h"
#include "pcg_basic.h"
#include <omp.h>

// Il rallentamento diffusivo scala come N^2. Z=2 è perfetto.
#define THERM_MIN_METROPOLIS 10000
#define THERM_MULT_METROPOLIS 20.0
#define THERM_Z_METROPOLIS 2

static int metropolis_therm_sweeps(int N) {
    double scaled = THERM_MULT_METROPOLIS * pow((double)N, THERM_Z_METROPOLIS);
    double n = scaled > THERM_MIN_METROPOLIS ? scaled : THERM_MIN_METROPOLIS;
    return (int)n;
}

// ----------------------------------------------------------------------------

void simulation(Polymer* poly, double eta, double beta, const char* alg, const int n_measures, const char* sim_name, pcg32_random_t* rng) {

    FILE* file = open_measurement_file(sim_name, poly->N, eta, beta);

    // The bar's row is the OpenMP thread slot this job landed on
    #ifdef _OPENMP
        int row = omp_get_thread_num();
    #else
        int row = 0;
    #endif

    char label[64];
    snprintf(label, sizeof(label), "N=%-4d eta=%.4f beta=%.4f", poly->N, eta, beta);

    if(strcmp(alg, "metropolis") == 0) {
        int N_therm = metropolis_therm_sweeps(poly->N);
        double max_delta = (0.5>2*sqrt(eta)) ? 2*sqrt(eta) : 0.5;

        // FASE 1: TERMALIZZAZIONE
        double therm_t0 = progress_now();
        int therm_step = N_therm / 100 > 0 ? N_therm / 100 : 1;
        
        for(int i = 0; i < N_therm; i++) {
            // 1. Mossa locale (diffusione di Trotter)
            metropolis_sweep(poly, eta, max_delta, rng, NULL);
            // 2. Mossa globale (campionamento dei settori topologici)
            //winding_move(poly, a, rng);

            // Controlla se i è multiplo di N_therm/100, o se ha finito
            if(i % therm_step == 0 || i == N_therm - 1) {
                progress_update(row, label, "therm", i + 1, N_therm, progress_now() - therm_t0);
            }
        }

        // FASE 2: MISURA

        poly->n_wrap = 0;
        poly->sum_d0 = 0;
        poly->sum_d1 = 0;

        double meas_t0 = progress_now();
        int meas_step = n_measures / 100 > 0 ? n_measures / 100 : 1;
        
        int acc_count = 0;

        for(int i = 0; i < n_measures; i++) {
            // Evoluzione del sistema
            metropolis_sweep(poly, eta, max_delta, rng, &acc_count);
            //winding_move(poly, a, rng);
            
            // Calcolo delle osservabili (l'energia dipende esplicitamente da 'a')
            double E = get_energy(poly, eta);
            int Q = get_winding(poly);
            
            // Salvataggio su file
            fprintf(file, "%d,%.10f,%d\n", i, E, Q);
            
            if(i % meas_step == 0 || i == n_measures - 1) {
                progress_update(row, label, "meas", i + 1, n_measures, progress_now() - meas_t0);
            }
        }

        FILE* file3 = fopen("wrap.txt", "a");
        fprintf(file3, "N=%d eta=%.4f  d*=%.4f  (dopo: %.4f)  n_wrap=%ld\n", poly->N, eta, poly->sum_d0/poly->n_wrap, poly->sum_d1/poly->n_wrap, poly->n_wrap);
        fclose(file3);

        double acc_rate = (double)acc_count/(double)(n_measures*poly->N);

        FILE* file2 = fopen("past.txt", "a");
        fprintf(file2, "eta: %.4f, max_delta: %.4f, acc_rate: %.4f\n", eta, max_delta, acc_rate*100);
        fclose(file2);

    } else {
        // Se in futuro vorrai implementare algoritmi avanzati per i Path Integral 
        // (es. Worm Algorithm o Staging), potrai aggiungerli qui.
        fprintf(stderr, "Errore: algoritmo %s non implementato per la particella su anello.\n", alg);
    }

    if (file != NULL) {
        fclose(file);
    }

    progress_job_done();
    return;
}


// ----------------------------------------------------------------------------
// Parallel Tempering

// Ogni quanti sweep proporre lo scambio
#define SWEEPS_PER_EXCHANGE 10

void simulation_pt(Polymer** replicas, double* eta_array, double* beta_array, int NR, int n_measures, const char* sim_name, pcg32_random_t* rngs) {
    FILE** files = malloc(NR * sizeof(FILE*));
    double max_delta = 0.5;
    
    for (int k = 0; k < NR; k++) {
        files[k] = open_measurement_file(sim_name, replicas[0]->N, eta_array[k], beta_array[k]);
    }

    int N_therm = metropolis_therm_sweeps(replicas[0]->N);
    
    // --- SETUP DELLA BARRA DI PROGRESSO ---
    int row = 0; // Il PT usa una singola riga unificata per tutte le repliche
    char label[64];
    snprintf(label, sizeof(label), "PT N=%-4d M=%-2d", replicas[0]->N, NR);
    
    double therm_t0 = progress_now();
    int therm_step = N_therm / 100 > 0 ? N_therm / 100 : 1;

    // FASE 1: TERMALIZZAZIONE
    for (int i = 0; i < N_therm; i++) {
        // L'opzione schedule (static) indica al compilatore di dividere il lavoro 
        // in blocchi rigidi e prestabiliti; essendo ogni polimero lungo esattamente N, 
        // il carico di lavoro è identico per tutti, rendendo questa la scelta più efficiente.
        #pragma omp parallel for schedule(static)
        for (int k = 0; k < NR; k++) {
            for (int s = 0; s < SWEEPS_PER_EXCHANGE; s++) {
                metropolis_sweep(replicas[k], eta_array[k], max_delta, &rngs[k], NULL);
            }
        }
        
        // Estrae 0 o 1 in modo casuale ad ogni ciclo di scambio, garantisce bilancio dettagliato ad ogni passo
        int offset = pcg32_random_r(&rngs[0]) % 2;

        for (int k = offset; k < NR - 1; k += 2) {
            replica_exchange(replicas[k], eta_array[k], replicas[k+1], eta_array[k+1], &rngs[0]);
        }

        // Aggiornamento grafico della termalizzazione
        if (i % therm_step == 0 || i == N_therm - 1) {
            progress_update(row, label, "therm", i + 1, N_therm, progress_now() - therm_t0);
        }
    }

    double meas_t0 = progress_now();
    int meas_step = n_measures / 100 > 0 ? n_measures / 100 : 1;

    for (int k = 0; k < NR; k++) {
        replicas[k]->sum_d0 = 0;
        replicas[k]->sum_d1 = 0;
    }

    // FASE 2: MISURA
    for (int i = 0; i < n_measures; i++) {
        #pragma omp parallel for schedule(static)
        for (int k = 0; k < NR; k++) {
            for (int s = 0; s < SWEEPS_PER_EXCHANGE; s++) {
                metropolis_sweep(replicas[k], eta_array[k], max_delta, &rngs[k], NULL);
            }
        }

        // Come prima: estrae 0 o 1 in modo casuale ad ogni ciclo di scambio, garantisce bilancio dettagliato ad ogni passo
        int offset = pcg32_random_r(&rngs[0]) % 2;

        for (int k = offset; k < NR - 1; k += 2) {
            replica_exchange(replicas[k], eta_array[k], replicas[k+1], eta_array[k+1], &rngs[0]);
        }

        #pragma omp parallel for schedule(static)
        for (int k = 0; k < NR; k++) {
            double E = get_energy(replicas[k], eta_array[k]);
            double Q = get_winding(replicas[k]);
            fprintf(files[k], "%d,%.10f,%.1f\n", i, E, Q);
            fflush(files[k]);
        }

        // Aggiornamento grafico delle misure
        if (i % meas_step == 0 || i == n_measures - 1) {
            progress_update(row, label, "meas", i + 1, n_measures, progress_now() - meas_t0);
        }
    }

    FILE* temp = fopen("P_acc.txt","a");
    fprintf(temp, "Nt: %d, C: %.4f, P_acc: %.4f\n", replicas[0]->N, (eta_array[1]/eta_array[0]), (double) replicas[0]->sum_d1/replicas[0]->sum_d0);
    fclose(temp);

    for (int k = 0; k < NR; k++) fclose(files[k]);
    free(files); 
}