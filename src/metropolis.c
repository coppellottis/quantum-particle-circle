#include "metropolis.h"
#include "pcg_basic.h"
#include <time.h>
#include <stdbool.h>
#include <math.h>
#include <stdio.h>

// Assumiamo massa m = 1, R=1 (Momento di inerzia I = 1)
void metropolis_update(Polymer* poly, double a, double max_delta, pcg32_random_t* rng, int i) {
    int N = poly->N;

    // Condizioni al contorno periodiche (anello chiuso)
    int left  = (i - 1 + N) % N; // +N evita segmentation fault
    int right = (i + 1) % N;

    double x_old = poly->beads[i];

    // 1. Proposta di spostamento continuo uniforme in [-max_delta, max_delta)
    double r_disp = (double)pcg32_random_r(rng) / 4294967296.0;
    double dx = max_delta * (2.0 * r_disp - 1.0); // proposta simmetrica rispetto a x_old
    double x_new = x_old + dx;

    // 2. Confinamento in (-pi, pi]
    if (x_new > 0.5)  x_new -= 1;
    if (x_new <= -0.5) x_new += 1;

    // 3. Calcolo della variazione dell'Azione (Delta S)
    // S_kin = (distanza)^2 / (2 * a)
    double d_left_old  = dist_circle(x_old, poly->beads[left]);
    double d_right_old = dist_circle(poly->beads[right], x_old);
    // Il fattore 4*pi^2 deriva dal riscalamento di x=phi/2pi
    double S_old = (4*M_PI*M_PI)*(d_left_old * d_left_old + d_right_old * d_right_old) / (2.0 * a);

    double d_left_new  = dist_circle(x_new, poly->beads[left]);
    double d_right_new = dist_circle(poly->beads[right], x_new);
    double S_new =  (4*M_PI*M_PI)*(d_left_new * d_left_new + d_right_new * d_right_new) / (2.0 * a);

    double delta_S = S_new - S_old;

    // In presenza di un potenziale:
    // delta_S += a * (V(x_new) - V(x_old));
    // Chiaramente, l'argomento del potenziale deve tenere conto del riscalamento /2pi

    // 4. Accettazione Metropolis (con cortocircuito per delta_S < 0)
    if (delta_S <= 0.0) {
        poly->beads[i] = x_new;
    } else {
        double t = (double)pcg32_random_r(rng) / 4294967296.0;
        // Chiamata costosa a exp() solo se strettamente necessario
        if (t < exp(-delta_S)) {
            poly->beads[i] = x_new;
        }
    }
}

void metropolis_sweep(Polymer* poly, double a, double max_delta, pcg32_random_t* rng) {
    int N = poly->N;
    for(int i=0; i < N; i++) {
        metropolis_update(poly, a, max_delta, rng, i);
    }
}

// Calcolo l'azione totale
static double get_action_sum(Polymer* poly) {
    double sum_du_sq = 0.0;
    for (int i = 0; i < poly->N; i++) {
        int next = (i + 1) % poly->N;
        double du = dist_circle(poly->beads[next], poly->beads[i]);
        sum_du_sq += (du * du);
    }
    return 4.0 * M_PI * M_PI * sum_du_sq;
}

void replica_exchange(Polymer* p1, double a1, Polymer* p2, double a2, pcg32_random_t* rng) {
    double K1 = get_action_sum(p1);
    double K2 = get_action_sum(p2);
    
    // Delta S = S(X1,a1)+S(X2,a2)-S(X2,a1)-S(X1,a1), semplificando
    // si ottiene l'espressione seguente
    double delta_S = (K2 - K1) * (1.0 / (2.0 * a1) - 1.0 / (2.0 * a2));
    
    // Scambio dei puntatori
    if (delta_S <= 0.0 || ((double)pcg32_random_r(rng) / 4294967296.0) < exp(-delta_S)) {
        double* temp_beads = p1->beads;
        p1->beads = p2->beads;
        p2->beads = temp_beads;
    }
}