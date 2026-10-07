#ifndef POLYMER_H
    #define POLYMER_H 1

    #include <stdbool.h>
    #include <stdint.h>
    #include <math.h>
    #include "pcg_basic.h"

    typedef struct{
        double* beads;
        int N;
    } Polymer;

    Polymer* init_polymer(int N, bool hot, pcg32_random_t* rng);
    void free_polymer(Polymer* poly);

    // Calcola la distanza geodetica più breve su un anello di circonferenza 1.0
    static inline double dist_circle(double u1, double u2) {
        double du = u1 - u2;
        // L'arrotondamento all'intero più vicino riavvolge automaticamente 
        // la distanza se supera mezzo giro (0.5)
        return du - round(du); 
    }

#endif
