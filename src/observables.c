#include <stdlib.h>
#include <math.h>
#include "observables.h"
#include "polymer.h"


// In order to avoid double counting, the sum
// over nearest neighbours is performed by considering only
// the right and the down neighbours. J=1.
double get_energy(Polymer* poly, double a) {
    int N = poly->N;
    double sum_du_sq = 0.0;
    
    for (int i = 0; i < N; i++) {
        int next = (i + 1) % N;
        double du = dist_circle(poly->beads[next], poly->beads[i]);
        sum_du_sq += (du * du);
    }
    
    double term_const = 1.0 / (2.0 * a);
    double term_fluct = (4.0 * M_PI * M_PI * sum_du_sq) / (2.0 * N * a * a);
    
    return (term_const - term_fluct);
}

double get_winding(Polymer* poly) {
    int N = poly->N;
    double winding_sum = 0.0;

    for (int i = 0; i < N; i++) {
        int next = (i + 1) % N;
        
        // Calcoliamo la distanza orientata (con segno) dal bead i al bead i+1.
        // È cruciale che dist_circle mantenga il segno (es. se passo da 0.4 a -0.4, 
        // la distanza più breve calcolata da dist_circle deve essere +0.2).
        double du = dist_circle(poly->beads[next], poly->beads[i]);
        
        winding_sum += du;
    }

    // Matematicamente winding_sum è già un intero esatto (es. 0.0, 1.0, -2.0).
    // Usiamo round() per correggere eventuali e minuscoli errori di 
    // arrotondamento accumulati dalla virgola mobile (es. 0.99999999999997).
    return round(winding_sum);
}