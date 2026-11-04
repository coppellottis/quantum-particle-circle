import numpy as np
import emcee
import pandas as pd
import  matplotlib.pyplot as plt

plt.rcParams.update({
    "font.family": "serif",
    "mathtext.fontset": "cm",
    "font.size": 20,

    "axes.labelsize": 20,
    "xtick.labelsize": 18,
    "ytick.labelsize": 18,

    "axes.linewidth": 1.0,

    "xtick.direction": "in",
    "ytick.direction": "in",

    "xtick.top": True,
    "ytick.right": True,

    "xtick.major.size": 6,
    "ytick.major.size": 6,
    "xtick.minor.size": 3,
    "ytick.minor.size": 3,

    "legend.fontsize": 20,
})

def block_obs(data, b_size, function, iterations = 100):
    b_number = len(data) // b_size
    blocks = data[:b_size*b_number].reshape(b_number, b_size)

    rng = np.random.default_rng()
    idx = rng.integers(0, b_number, size=(iterations, b_number))

    samples = blocks[idx].reshape(iterations, -1)
    estimates = np.array([function(sample) for sample in samples])

    return np.mean(estimates), np.std(estimates, ddof=1)

data = pd.read_csv("data/top_freezing/N300_a0.0067_beta2.0000.csv")
Q = np.array(data["Q"])

tau = emcee.autocorr.integrated_time(Q,2,1)[0]
print(len(Q)/tau)

kk = np.arange(1,2)*int(tau)

print(kk)

yy = []

for k in kk:
    Q2, Q2_std = block_obs(Q, k, lambda Q: Q**2)
    yy.append([Q2, Q2_std])

yy = np.array(yy)

fig, ax = plt.subplots(figsize=(12,8))

ax.plot(kk, yy[:,0],'-o')

ax2 = ax.twinx()
ax2.plot(kk, yy[:,1],'-o')

ax.set_xlabel(fr"Dimensione blocchi")
ax2.set_ylabel(fr"Errore stimato")
ax.set_ylabel(fr"$\langle Q^2\rangle$")

fig.tight_layout()

plt.show()