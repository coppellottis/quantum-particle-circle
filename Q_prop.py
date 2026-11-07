import numpy as np
import emcee
import pandas as pd
import random
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

def block_obs(data, b_size, function, iterations = 300):
    b_number = len(data) // b_size
    blocks = data[:b_size*b_number].reshape(b_number, b_size)

    rng = np.random.default_rng()

    # Un ricampionamento alla volta: in memoria c'e' sempre un solo campione
    # (len(data) elementi) invece di tutti i `iterations` insieme (iterations x len(data)).
    estimates = np.empty(iterations)
    for i in range(iterations):
        idx = rng.integers(0, b_number, size=b_number)
        sample = blocks[idx].reshape(-1)
        estimates[i] = np.mean(function(sample))   # function puo' restituire un array (es. Q**2): ne prendo la media

    return np.mean(estimates), np.std(estimates, ddof=1)

def b2(Q):
    Q = np.asarray(Q, dtype=float)          # evita problemi di dtype con Q**4
    q2 = np.mean(Q**2)
    q4 = np.mean(Q**4)
    return -(q4 - 3*q2**2) / (12*q2)

fig, ax = plt.subplots(figsize=(12,8))
N = [10]#np.arange(12,40,3)*10

yy_PT = []
yy_MT = []

for n in N:
    beta = 2
    eta = beta/n

    data = pd.read_csv(f"data/MT_0/N{n}_a{eta:.4f}_beta{beta:.4f}.csv", usecols=["Q"])
    Q = np.array(data["Q"])
    del data

    # mean = Q.mean()
    tau = emcee.autocorr.integrated_time(Q,10,1)[0]

    if np.isnan(tau):    # Q costante (congelata): tau non e' definito -> salto questa taglia
        print(f"PT  N{n}: tau=nan (Q costante?) -> salto\n")
        continue

    #_, std = block_obs(Q, int(tau), b2, iterations=100)
    #mean = b2(Q) 
    mean, std = block_obs(Q, int(tau), lambda Q: Q**2, iterations = 100)

    # std = Q.std()*np.sqrt(2*tau/len(Q))

    yy_PT.append([n, mean, std])

    print(f"N: {n}, <Q^2>: {mean}+-{std}\n")

    print(f"N: {n}, <Q>: {Q.mean()}+-{Q.std()*np.sqrt(2*tau/len(Q))}")
# for n in N:
#     beta = 2
#     eta = beta/n

#     data = pd.read_csv(f"data/met_bhat2/N{n}_a{eta:.4f}_beta{beta:.4f}.csv", usecols=["Q"])
#     Q = np.array(data["Q"])
#     del data

#     #mean = Q.mean()
#     tau = emcee.autocorr.integrated_time(Q,2,1)[0]

#     if np.isnan(tau):    # Q costante (congelata): tau non e' definito -> salto questa taglia
#         print(f"MT  N{n}: tau=nan (Q costante?) -> salto\n")
#         continue

#     _, std = block_obs(Q, int(tau), b2, iterations=100)
#     mean = b2(Q) 
#     #mean, std = block_obs(Q, int(tau), lambda Q: Q**2, iterations = 100)

#     # std = Q.std()*np.sqrt(2*tau/len(Q))

#     yy_MT.append([n, mean, std])

#     print(f"N{n}: n, <b2>: {mean}+-{std}\n")

# yy_PT = np.asarray(yy_PT)
# yy_MT = np.asarray(yy_MT)

# ax.errorbar(
#     yy_PT[:, 0],      # x
#     yy_PT[:, 1],      # y
#     yerr=yy_PT[:, 2], # errore su y

#     fmt='o',          # formato di marker/linea
#     # color='blue',     # colore principale
#     markersize=6,     # dimensione marker
#     # markerfacecolor='white',
#     # markeredgecolor='blue',

#     # ecolor='red',     # colore barre di errore
#     elinewidth=1.5,   # spessore barre
#     capsize=4,        # lunghezza dei "cappucci"
#     capthick=1.5,      # spessore cappucci

#     label="Parallel Tempering",
#     alpha=.5
# )

# ax.errorbar(
#     yy_MT[:, 0],      # x
#     yy_MT[:, 1],      # y
#     yerr=yy_MT[:, 2], # errore su y

#     fmt='o',          # formato di marker/linea
#     # color='blue',     # colore principale
#     markersize=6,     # dimensione marker
#     # markerfacecolor='white',
#     # markeredgecolor='green',

#     # ecolor='red',     # colore barre di errore
#     elinewidth=1.5,   # spessore barre
#     capsize=4,        # lunghezza dei "cappucci"
#     capthick=1.5,      # spessore cappucci

#     label="Metropolis",
#     alpha=.5
# )

# ax.set_xlabel(fr"$N_t$")
# ax.set_ylabel(fr"$b_2$")

# ax.axhline(y=0, color='r', linestyle='-.', label=fr"$b_2=0$")

# ax.legend(loc="best")
# ax.grid(alpha=.2)

# fig.tight_layout()

#fig.savefig("results/Q_b2_MTPT.pdf")
#plt.show()