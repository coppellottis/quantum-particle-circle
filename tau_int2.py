import numpy as np
import emcee
import pandas as pd
import  matplotlib.pyplot as plt
from scipy.optimize import curve_fit

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


def lin(x, m, q):
    return m*x+q

fig, ax = plt.subplots(figsize=(12,8))
N = np.arange(3,36,3)*10

yy = []

for n in N:
    beta = 2
    eta = beta/n

    data = pd.read_csv(f"data/tau_intb/N{n}_a{eta:.4f}_beta{beta:.4f}.csv")
    Q = np.array(data["Q"])
    rho = emcee.autocorr.function_1d(Q)

    tau, M = tau_int(rho, 2, n)

    yy.append([n, tau, M])

yy = np.array(yy)

ax.plot(yy[:,0], np.log(yy[:,1]), '-o',label=r"$\tau_{int}^Q (\Delta=0.5)$")

xx = np.linspace(min(yy[:,0]),max(yy[:,0]))

popt, pcov = curve_fit(lin, yy[:,0], np.log(yy[:,1]))

m, q = popt
dm, dq = np.sqrt(np.diag(pcov))

ax.plot(xx, lin(xx, m, q), '-.', color=ax.lines[-1].get_color(), label=fr"Fit: $e^{{aN_t}}, a\simeq {m:.3f}\pm{dm:.3f}$")

yy = []

for n in N:
    beta = 2
    eta = beta/n

    data = pd.read_csv(f"data/tau_int/N{n}_a{eta:.4f}_beta{beta:.4f}.csv")
    Q = np.array(data["Q"])
    rho = emcee.autocorr.function_1d(Q)

    tau, M = tau_int(rho, 2, n)

    yy.append([n, tau, M])

yy = np.array(yy)

ax.plot(yy[:,0], np.log(yy[:,1]), '-o',label=r"$\tau_{int}^Q (\Delta=\min(0.5,2\sqrt{\eta}))$")

xx = np.linspace(min(yy[:,0]),max(yy[:,0]))

popt, pcov = curve_fit(lin, yy[:,0], np.log(yy[:,1]))

m, q = popt
dm, dq = np.sqrt(np.diag(pcov))

ax.plot(xx, lin(xx, m, q), '-.', color=ax.lines[-1].get_color(), label=fr"Fit: $e^{{aN_t}}, a\simeq {m:.3f}\pm{dm:.3f}$")

ax.plot(xx, xx/16, 'r-.', label=r"$\propto e^{N_t/16}$")

ax.set_xlabel(fr"$N_t$")
ax.set_ylabel(fr"$\tau_{{int}}^{{Q}} [\log]$")

ax.legend(loc="best")

ax.grid(alpha=.2)

fig.tight_layout()
fig.savefig("results/tau_int2b.pdf", bbox_inches="tight")

plt.show()