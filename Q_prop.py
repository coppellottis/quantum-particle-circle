import numpy as np
import emcee
import pandas as pd
import random

def block_obs(data, b_size, function, iterations = 100):
    b_number = len(data) // b_size
    blocks = data[:b_size*b_number].reshape(b_number, b_size)

    rng = np.random.default_rng()
    idx = rng.integers(0, b_number, size=(iterations, b_number))

    samples = blocks[idx].reshape(iterations, -1)
    estimates = np.array([function(sample) for sample in samples])

    return np.mean(estimates), np.std(estimates, ddof=1)



data = pd.read_csv("data/var_10/N12_a0.1667_beta2.0000.csv")
Q = np.array(data["Q"])

mean = Q.mean()
tau = emcee.autocorr.integrated_time(Q,2,1)[0]
std = Q.std()*np.sqrt(2*tau/len(Q))

print(f"<Q>: {mean}+-{std}\n")

Q2, Q2_std = block_obs(Q, int(tau), lambda Q: Q**2)

print(f"<Q^2>: {Q2}+-{Q2_std}")