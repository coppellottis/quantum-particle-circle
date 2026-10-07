#ifndef SIMULATION_H
    #define SIMULATION_H 1

    #include "polymer.h"
    #include "pcg_basic.h"

    void simulation(Polymer* poly, double a, double beta, const char* alg, const int n_measures, const char* sim_name, pcg32_random_t* rng);
    void simulation_pt(Polymer** replicas, double* a_array, double* beta_array, int M, int n_measures, const char* sim_name, pcg32_random_t* rngs);

#endif