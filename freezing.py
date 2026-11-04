import numpy as np
import pandas as pd
import  matplotlib.pyplot as plt
import emcee

N = np.arange(1,5)*100
beta = 2

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

fig, ax = plt.subplots(4,1,figsize=(8,12),sharex=True)

for i, n in enumerate(N):
    eta = beta/n
    data = pd.read_csv(f"data/top_freezing/N{n}_a{eta:.4f}_beta{beta:.4f}.csv")

    ax.flat[i % 4].plot(data["sweep"], data["Q"],alpha=.4,label=rf"$N_t={n}$")

    mean = data["Q"].mean()
    tau = emcee.autocorr.integrated_time(data["Q"],2,1)
    err = data["Q"].std()*np.sqrt(2*tau/len(data["Q"]))
    
    ax.flat[i % 4].axhline(y=mean, color='red')
    ax.flat[i % 4].axhline(y=mean+err, color ='red', linestyle="-.", alpha=.5)
    ax.flat[i % 4].axhline(y=mean-err, color ='red', linestyle="-.", alpha=.5)
    #ax.flat[i % 4].set_xlabel("[sweep]")
    ax.flat[i % 4].set_ylabel("Q")
    ax.flat[i % 4].legend(loc="best")

ax[-1].set_xlabel("[sweep]")

fig.tight_layout()

fig.savefig(r"results/freezing.pdf",bbox_inches="tight")
plt.show()
