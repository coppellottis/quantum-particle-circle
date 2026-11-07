import numpy as np
import emcee
import pandas as pd
import  matplotlib.pyplot as plt

def tau_int(rho,c,n):
    tau = .5+np.cumsum(rho[1:])
    M = np.arange(1, len(rho))

    idx = np.where(M >= c*tau)[0]
    if len(idx) == 0:
        print(f"Nessuna finestra autoconsistente trovata (N_t:{n})")
        return 0,0

    M_star = M[idx[0]]
    tau_star = tau[idx[0]]

    return tau_star, M_star

# data = pd.read_csv("data/tau_int/N240_a0.0083_beta2.0000.csv")
# Q = np.array(data["Q"])

# rho = emcee.autocorr.function_1d(Q)
# tau, M = tau_int(rho, 2, 240)

# print(fr"tau: {tau}, N/tau={(1e6/tau)}")

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


fig, ax = plt.subplots(figsize=(12,8))
N = np.arange(12,40,3)*10


for n in N:
    beta = 2
    eta = beta/n

    data = pd.read_csv(f"data/PT_BHAT/N{n}_a{eta:.4f}_beta{beta:.4f}.csv")

    rho = emcee.autocorr.function_1d(np.array(data["Q"]))
    tau = 0.5 + np.cumsum(rho[1:])
    M = np.arange(1, len(rho))
    tau_star, M_star = tau_int(rho,10,n)


    ax.plot(M, tau,'-.',label=fr"$N_t={n}$")

    ax.scatter(M_star, tau_star,color=ax.lines[-1].get_color(), s=80, zorder=5)

###########################################

ax.set_xlabel(fr"$M$")
ax.set_ylabel(r"$\tau_{int}^{Q}(M)/\tau_{int}^{Q}(M^\star)$")

ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))

ax.grid(alpha=0.3)

ax.legend(loc="best")
ax.set_yscale("log")
ax.set_xscale("log")
ax.set_ylim(10**(-.2),1e2)
ax.set_xlim(1e0,1e3)

fig.tight_layout()


fig.savefig("results/tau_int_PT.pdf")
plt.show()