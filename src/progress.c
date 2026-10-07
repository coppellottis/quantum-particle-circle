#include <stdio.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include "progress.h"

static int g_n_rows = 0;    // per-job bar rows (one per worker thread)
static int g_n_jobs = 0;    // total jobs in the batch
static int g_jobs_done = 0;
static double g_batch_start = 0.0;
static int g_is_tty = 0;

double progress_now(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9;
}

// Formats a duration in seconds as e.g. "42s", "3m07s", "1h23m". A
// negative value (unknown/not-yet-estimable) prints as "--".
static void format_duration(double seconds, char* buf, int bufsize) {
    if (seconds < 0) {
        snprintf(buf, bufsize, "--");
        return;
    }
    long s = (long)(seconds + 0.5);
    if (s < 60) {
        snprintf(buf, bufsize, "%lds", s);
    } else if (s < 3600) {
        snprintf(buf, bufsize, "%ldm%02lds", s / 60, s % 60);
    } else {
        snprintf(buf, bufsize, "%ldh%02ldm", s / 3600, (s % 3600) / 60);
    }
}

static void draw_bar(char* buf, int bufsize, double frac) {
    const int width = 20;
    if (frac < 0) frac = 0;
    if (frac > 1) frac = 1;
    int filled = (int)(frac * width);

    int pos = 0;
    buf[pos++] = '[';
    for (int i = 0; i < width && pos < bufsize - 1; i++) {
        buf[pos++] = (i < filled) ? '#' : '-';
    }
    if (pos < bufsize - 1) buf[pos++] = ']';
    buf[pos] = '\0';

    char pct[16];
    snprintf(pct, sizeof(pct), " %5.1f%%", frac * 100.0);
    strncat(buf, pct, bufsize - strlen(buf) - 1);
}

// Must be called from inside the "sim_progress" critical section
// (both progress_init and progress_job_done draw the header).
static void redraw_header_locked(int done) {
    double elapsed = progress_now() - g_batch_start;
    char elapsed_str[64], eta_str[64], header[256];

    format_duration(elapsed, elapsed_str, sizeof(elapsed_str));

    if (done > 0 && done < g_n_jobs) {
        // Rough estimate: average wall-clock time per completed job so
        // far, scaled down by how many jobs run concurrently. Real job
        // cost varies with L and beta (a big-L job near beta_c is much
        // more expensive than a small-L one far from criticality), so
        // this drifts as the batch works through jobs of very different
        // sizes -- but it's the only data available without knowing every
        // future job's cost up front, and it's enough to answer "roughly
        // how much longer".
        double avg_per_job = elapsed / done;
        double eta = avg_per_job * (g_n_jobs - done) / g_n_rows;
        format_duration(eta, eta_str, sizeof(eta_str));
        snprintf(header, sizeof(header), "Job completati: %d/%d | trascorsi: %s | stima rimanente: ~%s",
                 done, g_n_jobs, elapsed_str, eta_str);
    } else if (done >= g_n_jobs && g_n_jobs > 0) {
        snprintf(header, sizeof(header), "Job completati: %d/%d | tempo totale: %s", done, g_n_jobs, elapsed_str);
    } else {
        snprintf(header, sizeof(header), "Job completati: %d/%d | trascorsi: %s | stima rimanente: --",
                 done, g_n_jobs, elapsed_str);
    }

    if (g_is_tty) {
        int up = g_n_rows + 1; // header sits above all g_n_rows bar rows
        printf("\033[%dA\r\033[K%s", up, header);
        printf("\033[%dB\r", up);
        fflush(stdout);
    } else {
        printf("%s\n", header);
    }
}

void progress_init(int n_rows, int n_jobs) {
    g_n_rows = (n_rows > 0) ? n_rows : 1;
    g_n_jobs = (n_jobs > 0) ? n_jobs : 1;
    g_jobs_done = 0;
    g_batch_start = progress_now();
    g_is_tty = isatty(fileno(stdout));

    if (g_is_tty) {
        // Reserve g_n_rows blank bar lines plus one header line above
        // them. Every later update moves the cursor up from this
        // baseline to the right row, redraws it, then moves back down, so
        // it never assumes an absolute screen position -- only a
        // relative offset from this baseline.
        for (int i = 0; i < g_n_rows + 1; i++) {
            printf("\n");
        }
        fflush(stdout);

        #pragma omp critical(sim_progress)
        {
            redraw_header_locked(0);
        }
    }
}

void progress_update(int row, const char* label, const char* phase, long done, long total, double elapsed_sec) {
    if (row < 0) row = 0;
    if (g_n_rows > 0 && row >= g_n_rows) row = g_n_rows - 1;

    double frac = (total > 0) ? (double)done / (double)total : 0.0;
    char bar[40];
    draw_bar(bar, sizeof(bar), frac);

    char eta_str[64];
    if (done >= total) {
        snprintf(eta_str, sizeof(eta_str), "0s");
    } else if (elapsed_sec > 0.05 && done > 0) {
        double rate = (double)done / elapsed_sec; // items/sec
        double eta = (double)(total - done) / rate;
        format_duration(eta, eta_str, sizeof(eta_str));
    } else {
        snprintf(eta_str, sizeof(eta_str), "--");
    }

    // The whole "move cursor, clear, redraw, move back" sequence has to
    // land on the terminal as one atomic unit, otherwise two threads
    // updating different rows at the same time can interleave their
    // escape codes and text into garbage. A critical section is cheap
    // here: updates are throttled by the caller to ~100 per phase, not
    // once per sweep.
    #pragma omp critical(sim_progress)
    {
        if (g_is_tty) {
            int up = g_n_rows - row;
            printf("\033[%dA\r\033[K%-30s %-6s %s eta %s", up, label, phase, bar, eta_str);
            printf("\033[%dB\r", up);
            fflush(stdout);
        } else {
            printf("%-30s %-6s %s eta %s\n", label, phase, bar, eta_str);
        }
    }
}

void progress_job_done(void) {
    #pragma omp critical(sim_progress)
    {
        g_jobs_done++;
        redraw_header_locked(g_jobs_done);
    }
}

void progress_finish(void) {
    if (g_is_tty) {
        printf("\n");
        fflush(stdout);
    }
}
