import inspect
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


rng = np.random.default_rng()

# Binder's cumulant
def binder_cumulant(x) :
    return 1-(x**4).mean()/(3*((x**2).mean())**2)

# block bootstrap error
# use this for secondary variables such as U4 Binder cumulant
# returns error for different block dimension k
def blckbstr_error(x, f, k = 20):
    k = max(1, int(k))

    N = x.shape[0]
    n = N // k

    # pre-computation of block
    blocks = np.array([x[i*k:(i+1)*k] for i in range(n)])

    R = 1000
    means = np.zeros(R)

    for r in range(R):
        indices = rng.integers(0, n, size=n)
        means[r] = f(blocks[indices].reshape(-1))

    error = means.std(ddof=1)

    return error

# original definition of the integrated autocorrelation time.
# O(N^2): slow, only meant as a reference/cross-check on small samples
# for tau_int_fft below. Uses the same Sokal automatic-windowing
# criterion (same default c) so the two are directly comparable --
# they used to disagree (c=10 here vs c=5 in tau_int_fft), which made
# them silently report different tau for the same data depending on
# which one you called.
def tau_int(x, c=5) :
    x = np.asarray(x)
    tau = 0.5
    N = x.shape[0]

    mean = x.mean()
    var = x.var()

    for k in range(1,N//2) :
        d = var*(N-k)
        rho_k = 0
        for i in range(N-k):
            rho_k += (x[i]-mean)*(x[i+k]-mean)

        tau += 1/d*rho_k

        if(k >= c*tau) :
            return tau

    print("Warning: Window did not converge; returned None.")
    return None

# Running integrated-autocorrelation-time estimate tau(M) for every window
# M = 0, 1, ..., N-1 (plus the underlying normalized autocorrelation
# function rho), computed via FFT. This is the same autocovariance
# tau_int_fft has always used internally -- pulled out into its own
# function so a diagnostic plot of tau(M) vs M (see tau_diagnostic.py) can
# be drawn from exactly the numbers tau_int_fft itself works with, instead
# of a second, separately-written computation that could silently drift
# out of sync with it.
def tau_running(x):
    x = np.asarray(x)
    N = len(x)

    x = x - np.mean(x)

    nfft = 2**int(np.ceil(np.log2(2*N)))

    f = np.fft.fft(x, n=nfft)
    power = f * np.conjugate(f)

    acov = np.fft.ifft(power).real[:N]
    acov /= N # keeps only the first 0,... N-1 elements (on nnft)
    # biased convention used; acov/= N works better than un-biased acov/= N-lag... (biased+more stable for large k)

    rho = acov / acov[0]

    taus = np.concatenate(([0.5], 0.5+np.cumsum(rho[1:])))
    return taus, rho

# Sokal's automatic windowing: the smallest window M for which
# M >= c*tau(M) (self-consistent stopping criterion), subject to a minimum
# window and to M staying below N/2 (a window that needs half the data to
# "converge" isn't trustworthy -- it means N is too small relative to tau,
# most often because of critical slowing down at large L close to beta_c,
# or insufficient thermalization). Returns the chosen M (an index into
# taus), or None if no M satisfies the criterion within N/2. Split out of
# tau_int_fft for the same reason as tau_running above: tau_diagnostic.py
# needs to know exactly which M was picked, not just the resulting tau.
def sokal_window(taus, N, c=5, min_window=4):
    idx = np.arange(0, len(taus)) # lag k
    valid = (idx >= min_window) & (idx < N // 2) & (idx >= c * taus)
    if not np.any(valid):
        return None
    return int(idx[np.argmax(valid)])

# (faster) computation of the int. autocorr. time through fast fourier transform.
#
# label: optional free-form string (e.g. "wolff L=64 beta=0.4407 obs=|m|")
# included in the non-convergence warning, so that when this is called in
# a loop over many (L, beta, observable) combinations -- as analyze.py and
# tau.py do -- you can actually tell *which* run needs more statistics or
# better thermalization instead of an unattributed "Warning: Window did
# not converge" printed somewhere in a wall of output.
def tau_int_fft(x, c = 6, min_window = 4, label = None) :

    x = np.asarray(x)
    N = len(x)

    taus, rho = tau_running(x)
    M = sokal_window(taus, N, c=c, min_window=min_window)

    if M is None:
        tag = f" [{label}]" if label else ""
        print(f"Warning: autocorrelation window did not converge{tag} "
              f"(N={N} samples not >> tau; more measurements or longer "
              f"thermalization needed). Returned None.")
        return None
    else :
        return max(taus[M], 0.5), M

def error(x, label=None) :
    x = np.asarray(x)
    N = len(x)

    tau, _ = tau_int_fft(x, label=label)
    if tau is None:
        return np.nan, tau
    else:
        error = np.sqrt(np.var(x, ddof=1)*2*tau/N)
        return error, tau

# Calls tau_func(x, label=label) if tau_func accepts a label kwarg (as
# tau_int_fft does, for a labeled non-convergence warning), otherwise
# falls back to tau_func(x). Keeps get_tau generic over any tau_func while
# still propagating diagnostics through to tau_int_fft when possible,
# instead of printing one unlabeled warning from inside tau_func and a
# second, differently-labeled one from get_tau itself. Checking the
# signature (rather than try/except TypeError around the call) avoids
# accidentally swallowing a genuine TypeError raised from inside tau_func.
def _call_tau_func(tau_func, x, label):
    try:
        accepts_label = 'label' in inspect.signature(tau_func).parameters
    except (TypeError, ValueError):
        accepts_label = False
    if accepts_label:
        return tau_func(x, label=label)
    return tau_func(x)

# tau value + block jackknife error
def get_tau(x, tau_func, n_blocks=20, label=None):
    x = np.asarray(x)
    N = len(x)

    tau, M = _call_tau_func(tau_func, x, label)
    if tau is None:
        tag = f" [{label}]" if label else ""
        print(f"Warning: get_tau could not get a point estimate for tau{tag}; returning (None, nan).")
        return None, np.nan

    tau_err = np.sqrt(2*(2*M+1-tau)/N*(tau**2))
    
    return tau, tau_err

# Safe block size (in samples) for blckbstr_error() from a possibly-missing
# tau estimate. analyze.py used to pass np.floor(10*tau_x) straight in;
# when tau_x was None (tau_int_fft didn't converge for that run) this
# raised a TypeError/ValueError and killed the whole analyze() job for
# that L, silently dropping every observable for that lattice size, not
# just the one that failed. Falling back to N/fallback_frac keeps going
# with a reasonable (if less well-motivated) block size instead.
def safe_block_k(tau, N, fallback_frac=20, label=None):
    if tau is None or not np.isfinite(tau):
        tag = f" [{label}]" if label else ""
        print(f"Warning: no valid tau for block-bootstrap{tag}; falling back to N/{fallback_frac} as block size.")
        return max(1, N // fallback_frac)
    return max(1, int(np.floor(10*tau)))

def build_histogram(energies, e_min, e_max, de=4) :
    # note that the sup limit is exclued while using np.arange
    # (reason why +4 must be added to np.max(energies))
    energies = np.asarray(energies, dtype=float)
    E_bins = np.arange(e_min, e_max + de, de)
    H = np.zeros(len(E_bins), dtype=int)
    idx = np.round((energies - e_min) / 4).astype(int)
    H = np.bincount(idx, minlength=len(E_bins))[:len(E_bins)]

    return E_bins, H

# output format
def format_error(value, error, sig=1):
    """
    Format as value(error), e.g.
    2.1734 ± 0.0512 -> 2.17(5)
    2.1734 ± 0.0051 -> 2.173(5)
    """
    if error == 0:
        return f"{value}"

    exponent = int(np.floor(np.log10(abs(error))))
    decimals = max(0, -exponent + (sig - 1))

    value_r = round(value, decimals)
    error_r = round(error, decimals)

    # errore espresso come intero nelle ultime cifre
    error_digits = int(round(error_r * 10**decimals))

    return f"{value_r:.{decimals}f}({error_digits})"
