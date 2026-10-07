#ifndef PROGRESS_H
    #define PROGRESS_H 1

    // Monotonic clock, in seconds (CLOCK_MONOTONIC -- immune to wall-clock
    // adjustments, unlike time()). Callers use this to time how long a
    // phase (thermalization/measurement) has been running, which is what
    // progress_update() needs to extrapolate an ETA.
    double progress_now(void);

    // Reserves n_rows per-job bar rows plus one header line above them
    // for an overall batch summary (jobs completed, elapsed, ETA), and
    // starts the batch clock used for that summary. n_rows is one slot
    // per OpenMP thread actually used to run jobs, not one per job: with
    // a large scan there can be far more (L, beta) jobs than CPU cores,
    // and a bar per queued job would print hundreds of mostly-idle lines.
    // n_jobs is the total number of jobs in the batch, used for the
    // header's "X/N completed" and ETA. Call once, from a single thread,
    // before entering the parallel region.
    void progress_init(int n_rows, int n_jobs);

    // Updates the bar in slot "row" (0..n_rows-1, typically
    // omp_get_thread_num()) to show "label" (the job's L/beta/algo), the
    // current "phase" ("therm"/"meas"), done/total progress within that
    // phase, and an ETA for the phase extrapolated from elapsed_sec (time
    // spent in the phase so far, from progress_now() at the phase's
    // start) and the observed rate done/elapsed_sec. Safe to call
    // concurrently from multiple threads (updates are internally
    // serialized). When stdout isn't a terminal (redirected to a file or
    // piped), falls back to plain, non-overlapping progress lines instead
    // of cursor-positioning escape codes.
    void progress_update(int row, const char* label, const char* phase, long done, long total, double elapsed_sec);

    // Marks one job as finished (thread-safe) and refreshes the header
    // line with jobs-done/total, elapsed batch time, and an estimated
    // time remaining for the *whole* batch -- extrapolated from the
    // average wall-clock time per completed job so far, divided across
    // the number of worker rows (a rough estimate: real job cost varies
    // with L and beta, but it's the only data available without knowing
    // every future job's cost in advance). Call once per finished job,
    // from any thread.
    void progress_job_done(void);

    // Moves the cursor past the reserved bar area, leaving the terminal
    // clean for whatever prints next. Call once, from a single thread,
    // after the parallel region ends.
    void progress_finish(void);
#endif
