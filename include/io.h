#ifndef IO_H
    #define IO_H 1

    #include <stdio.h>

    // Writes data/<sim_name>/metadata.csv describing every (L, beta) job
    // in a simulation batch. Exits the process with a clear error message
    // if the file can't be opened (the previous inline version in main.c
    // never checked fopen()'s return value, so a missing/unwritable
    // data/<sim_name> folder would crash later with an obscure
    // segfault from fprintf(NULL, ...) instead).
    void write_metadata(const char* sim_name, int phase, double fixed_val, 
                    const int* N_array, int n_N, int n_measures);

    // Opens data/<sim_name>/<alg>_L<L>_beta<beta>.csv for writing and
    // writes its header row. The caller owns the returned FILE* and must
    // fclose() it.
    FILE* open_measurement_file(const char* sim_name, int N, double a, double beta);

#endif
