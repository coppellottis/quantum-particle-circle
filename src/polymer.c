#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include "polymer.h"
#include "pcg_basic.h"


Polymer* init_polymer(int N, bool hot, pcg32_random_t* rng) {
    Polymer* poly = malloc(sizeof(Polymer));
    poly->beads = malloc(N*sizeof(double));
    poly->N = N;

    for (int i = 0; i < N; i++) {
        if (hot) {
            // Hot Start: Posizioni casuali uniformemente distribuite in [-pi, pi)
            // Dividiamo per 2^32 per ottenere un double in [0, 1) dal PCG32
            double rand_01 = (double)pcg32_random_r(rng) / 4294967296.0;
            poly->beads[i] = -0.5 + 2.0 * 0.5 * rand_01; // variabile adimensionale x=theta/pi
        } else {
            // Cold Start: Tutti i bead concentrati nell'origine
            poly->beads[i] = 0.0;
        }
    }
    return poly;
}


void free_polymer(Polymer* poly) {
    free (poly);
    return;
}
