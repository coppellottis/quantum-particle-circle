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
///_______________________

int main(void){

    int n_N = 0; // numero di taglie della catena
    double I = 1.0; // momento d'inerzia
    double theta = 0; // vuoto theta
    
    char alg[] = "metropolis";
    char sim_name[256];


    // IMPOSTAZIONI GENERALI SIMULAZIONE
    printf("Nome della simulazione: ");
    scanf("%255s", sim_name);

    char folder_path[512];

    snprintf(folder_path, sizeof(folder_path), "data/%s", sim_name);

    if (mkdir("data", 0777) != 0 && errno != EEXIST) {
        perror("WARNING: Non è stato possibile creare la cartella data/.");
        exit(EXIT_FAILURE);
    }

    if (mkdir(folder_path, 0777) != 0) {
        if(errno == EEXIST) {
            perror("WARNING: La cartella già esiste. Proseguendo è possibile che i precedenti file vengano rimossi o sostituiti.");
        } else {
            perror("WARNING: Non è stato possibile creare la cartella della simulazione.");
            exit(EXIT_FAILURE);
        }
    }

    // SCELTA DELLA FASE DI LAVORO
    int phase = 0;
    while(phase < 1 || phase > 3) {
        printf("\n=== SELEZIONE DEL FLUSSO DI LAVORO ===\n");
        printf("1. Fase 1: Ricerca del Ground State (Fissa 'eta', esplora 'N')\n");
        printf("2. Fase 2: Limite del Continuo (Fissa 'beta', esplora 'N')\n");
        printf("3. Fase 3: Parallel Tempering\n");
        printf("Scegli la fase (1, 2 o 3): ");
        scanf("%d", &phase);
    }

    double fixed_val = 0.0;

    if (phase == 1) {
        printf("\n--- FASE 1: RICERCA DEL GROUND STATE ---\n");
        while(fixed_val <= 0.0) {
            printf("Inserisci il passo temporale fisso 'eta' (es. 0.1): ");
            scanf("%lf", &fixed_val);
        }
    } else if(phase == 2){
        printf("\n--- FASE 2: LIMITE DEL CONTINUO ---\n");
        while(fixed_val <= 0.0) {
            printf("Inserisci il tempo immaginario fisso 'beta' (='N*eta') (es. 2.0): ");
            scanf("%lf", &fixed_val);
        }
    } else if (phase == 3) {
        int Nt, NR, n_measures;
        double C, eta0;

        printf("\n--- FASE 3: PARALLEL TEMPERING ---\n");
        printf("Inserisci N (numero di bead fisso): ");
        scanf("%d", &Nt);
        printf("Inserisci M (numero di repliche/temperature): ");
        scanf("%d", &NR);
        printf("Inserisci parametro C: ");
        scanf("%lf", &C);
        printf("Inserisci eta iniziale: ");
        scanf("%lf", &eta0);
        double input = 0.0;
        printf("Inserisci il numero di misure: ");
        scanf("%lf", &input);
        if (input > 0) {
            n_measures = (int) input;
        }

        Polymer** replicas = malloc(NR * sizeof(Polymer*));
        double* eta_array = malloc(NR * sizeof(double));
        double* beta_array = malloc(NR * sizeof(double));
        pcg32_random_t* rngs = malloc(NR * sizeof(pcg32_random_t));

        for (int k = 0; k < NR; k++) {
            eta_array[k] = pow(C,k)*eta0;
            beta_array[k] = eta_array[k]*Nt;
            pcg32_srandom_r(&rngs[k], 42u ^ (uint64_t)k, 54u ^ (uint64_t)Nt);
            bool hot_start = true;
            replicas[k] = init_polymer(Nt, hot_start, &rngs[k]);
        }

        // --- SALVATAGGIO METADATI PER ANALYZE.PY ---
        write_metadata_pt(sim_name, Nt, NR, eta_array, beta_array, n_measures);

        omp_set_num_threads(NR);

        printf("\nAvvio simulazione PT con %d repliche...\n", NR);
        progress_init(1, 1);

        simulation_pt(replicas, eta_array, beta_array, NR, n_measures, sim_name, rngs);

        progress_job_done();
        progress_finish();

        for(int k = 0; k < NR; k++) free_polymer(replicas[k]);
        free(replicas); 
        free(eta_array); 
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
        if(i>0) {
            if(phase == 1) printf("[eta=%3.f, N=%d] Valore di beta corrispondente (beta = N*eta): %3.f\n", fixed_val, N_array[i-1], (double)N_array[i-1]*fixed_val);
        }
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
        double input = 0.0;
        printf("Inserisci il numero di sweep di misura: ");
        scanf("%lf", &input);
        if (input > 0) {
            n_measures = (int) input;
            break;
        }
    }
    
    // Riepilogo a schermo
    printf("\n=== SETUP SIMULAZIONE COMPLETATO ===\n");
    if (phase == 1) {
        printf("Modalita': Ricerca Ground State (eta fisso = %3.f)\n", fixed_val);
    } else {
        printf("Modalita': Limite Continuo (beta fisso = %3.f)\n", fixed_val);
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
        double eta, beta;

        if (phase == 1) {
            eta = fixed_val;          
            beta = N * eta;           
        } else {
            beta = fixed_val;       
            eta = beta / (double)N;   
        }

        pcg32_random_t rng;
        pcg32_srandom_r(&rng, 42u ^ (uint64_t)i, 54u ^ (uint64_t)N);

        bool hot_start = false; 
        Polymer* poly = init_polymer(N, hot_start, &rng);

        simulation(poly, eta, beta, alg, n_measures, sim_name, &rng);

        free_polymer(poly);
    }

    // Chiude in modo pulito l'interfaccia grafica del terminale
    progress_finish(); 
    free(N_array);

    return 0;
}