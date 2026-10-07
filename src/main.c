#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <errno.h>
#include <omp.h>
#include "pcg_basic.h"
#include "polymer.h"
#include "metropolis.h"
#include "io.h"
#include "simulation.h"
#include "observables.h"
#include "progress.h"

int main(void){

    int n_N = 0; // numero di taglie della catena
    double I = 1.0; // momento d'inerzia
    double theta = 0; // vuoto theta
    
    char alg[] = "metropolis";
    char sim_name[256];

    printf("Nome della simulazione: ");
    scanf("%255s", sim_name);

    char folder_path[512];

    snprintf(folder_path, sizeof(folder_path), "data/%s", sim_name);

    if (mkdir("data", 0777) != 0 && errno != EEXIST) {
        perror("WARNING: Non è stato possibile creare la cartella data/.");
        exit(EXIT_FAILURE);
    }

    if (mkdir(folder_path, 0777) != 0) {
        perror("WARNING: Non è stato possibile creare la cartella della simulazione.");
        exit(EXIT_FAILURE);
    }

    // 1. Scelta della Fase del Flusso di Lavoro
    int phase = 0;
    while(phase < 1 || phase > 3) {
        printf("\n=== SELEZIONE DEL FLUSSO DI LAVORO ===\n");
        printf("1. Fase 1: Ricerca del Ground State (Fissa 'a', esplora 'N')\n");
        printf("2. Fase 2: Limite del Continuo (Fissa 'beta', esplora 'N')\n");
        printf("3. Fase 3: Parallel Tempering (Fissa 'N', esplora range di 'beta')\n");
        printf("Scegli la fase (1, 2 o 3): ");
        scanf("%d", &phase);
    }

    double fixed_val = 0.0;

    if (phase == 1) {
        printf("\n--- FASE 1: RICERCA DEL GROUND STATE ---\n");
        while(fixed_val <= 0.0) {
            printf("Inserisci il passo temporale fisso 'a' (es. 0.1): ");
            scanf("%lf", &fixed_val);
        }
    } else if(phase == 2){
        printf("\n--- FASE 2: LIMITE DEL CONTINUO ---\n");
        while(fixed_val <= 0.0) {
            printf("Inserisci il tempo immaginario fisso 'beta' (='N*a') (es. 5.0): ");
            scanf("%lf", &fixed_val);
        }
    } else if (phase == 3) {
        int N, M, n_measures;
        double beta_min, beta_max;

        printf("\n--- FASE 3: PARALLEL TEMPERING ---\n");
        printf("Inserisci N (numero di bead fisso): ");
        scanf("%d", &N);
        printf("Inserisci M (numero di repliche/temperature): ");
        scanf("%d", &M);
        // C'è un modo migliore per determinare i beta (per massimizzare le accettazioni degli scambi...)
        printf("Inserisci beta minimo e massimo (es. 0.5 15.0): ");
        scanf("%lf %lf", &beta_min, &beta_max);
        printf("Inserisci il numero di misure: ");
        scanf("%d", &n_measures);

        Polymer** replicas = malloc(M * sizeof(Polymer*));
        double* a_array = malloc(M * sizeof(double));
        double* beta_array = malloc(M * sizeof(double));
        pcg32_random_t* rngs = malloc(M * sizeof(pcg32_random_t));

        for (int k = 0; k < M; k++) {
            beta_array[k] = beta_min + k * (beta_max - beta_min) / (M > 1 ? M - 1 : 1);
            a_array[k] = beta_array[k] / (double)N;
            pcg32_srandom_r(&rngs[k], 42u ^ (uint64_t)k, 54u ^ (uint64_t)N);
            replicas[k] = init_polymer(N, false, &rngs[k]);
        }

        omp_set_num_threads(M);

        printf("\nAvvio simulazione PT con %d repliche...\n", M);
        simulation_pt(replicas, a_array, beta_array, M, n_measures, sim_name, rngs);

        for(int k = 0; k < M; k++) free_polymer(replicas[k]);
        free(replicas); 
        free(a_array); 
        free(beta_array); 
        free(rngs);
        
        return 0; // Termina qui per la Fase 3
    }


    // 2. Acquisizione dei valori di N (chiamato L_array per coerenza con il tuo codice)
    while(1) {
        printf("Quante taglie di N (numero di bead) vuoi simulare? ");
        scanf("%d", &n_N);
        if (n_N > 0) break;
        printf("ATTENZIONE: Il numero di taglie deve essere > 0.\n");
    }

    int *N_array = (int *)malloc(n_N * sizeof(int));
    
    for (int i = 0; i < n_N; i++) {
        while(1) {
            printf("Inserisci N per il reticolo %d: ", i + 1);
            scanf("%d", &N_array[i]);
            if (N_array[i] > 0) break;
            printf("ATTENZIONE: N deve essere > 0.\n");
        }
    }

    // 3. Acquisizione dei parametri Monte Carlo
    int n_measures;
    
    printf("\n--- PARAMETRI MONTE CARLO ---\n");

    while(1) {
        printf("Inserisci il numero di sweep di misura: ");
        scanf("%d", &n_measures);
        if (n_measures > 0) break;
    }

    // Riepilogo a schermo
    printf("\n=== SETUP SIMULAZIONE COMPLETATO ===\n");
    if (phase == 1) {
        printf("Modalita': Ricerca Ground State (a fisso = %f)\n", fixed_val);
    } else {
        printf("Modalita': Limite Continuo (beta fisso = %f)\n", fixed_val);
    }
    printf("Cartella di output: data/%s\n", sim_name);
    printf("Taglie (N): %d\n", n_N);
    printf("====================================\n\n");
    
    write_metadata(sim_name, phase, fixed_val, N_array, n_N, n_measures);

    // Acquisizione thread OpenMP
    int n_threads;
    while(1) {
        printf("Inserisci il numero di thread OpenMP da usare: ");
        scanf("%d", &n_threads);
        if (n_threads > 0) break;
    }
    omp_set_num_threads(n_threads);

    // Inizializzazione della barra di progresso (se usi progress.h)
    // progress_init(n_threads, n_N);

    // =========================================================================
    // ESECUZIONE PARALLELA DELLE SIMULAZIONI
    // =========================================================================
    
    // Inizializzazione della barra di progresso: 
    // Alloca n_threads righe sul terminale per stampare le barre simultaneamente
    progress_init(n_threads, n_N);

    // =========================================================================
    // ESECUZIONE PARALLELA DELLE SIMULAZIONI
    // =========================================================================
    
    #pragma omp parallel for schedule(dynamic, 1)
    for (int i = 0; i < n_N; i++) {
        int N = N_array[i];
        double a, beta;

        if (phase == 1) {
            a = fixed_val;          
            beta = N * a;           
        } else {
            beta = fixed_val;       
            a = beta / (double)N;   
        }

        pcg32_random_t rng;
        pcg32_srandom_r(&rng, 42u ^ (uint64_t)i, 54u ^ (uint64_t)N);

        bool hot_start = false; 
        Polymer* poly = init_polymer(N, hot_start, &rng);

        simulation(poly, a, beta, alg, n_measures, sim_name, &rng);

        free_polymer(poly);
    }

    // Chiude in modo pulito l'interfaccia grafica del terminale
    progress_finish(); 
    free(N_array);

    return 0;
}