#ifndef METROPOLIS_H
    #define METROPOLIS_H 1

    #include "polymer.h"
    #include "pcg_basic.h"

    void metropolis_update(Polymer* poly, double a, double max_delta, pcg32_random_t* rng, int i);
    void metropolis_sweep(Polymer* poly, double a, double max_delta, pcg32_random_t* rng);
    static double get_action_sum(Polymer* poly);
    void replica_exchange(Polymer* p1, double a1, Polymer* p2, double a2, pcg32_random_t* rng);
    
#endif