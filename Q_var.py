import numpy as np
import  matplotlib.pyplot as plt
import math

def Q_var(beta) :
    Q = np.arange(-10000,+10001)
    w = np.exp(-Q[:,None]**2/(2*beta[None,:]))
    den = np.sum(w,axis=0)
    num = np.sum((Q[:,None]**2)*w,axis=0)
    return num/den

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

fig, ax = plt.subplots(figsize=(12, 8))

betas = np.linspace(1e-1,5,1000)
var = Q_var(betas)

print(Q_var(np.array([1/(2*math.pi)])))
print((1/(2*math.pi))-Q_var(np.array([1/(2*math.pi)])))

ax.plot(betas,var,label=r"$\langle Q^2\rangle$",linewidth=3)

ax.plot(betas,2*np.exp(-1/(2*betas)),label=r"$\langle Q^2\rangle=2e^{-1/2\hat\beta}$",linestyle='--')

ax.axline((0, 0), slope=1, color='red', linestyle='--', label=r'$\langle Q^2\rangle = \hat{\beta}$')
ax.axvline(x=2,color='black', linestyle="-.",label=r"$\hat\beta=2$")

ax.set_xlabel(rf'$\hat\beta$')
ax.set_ylabel(rf"$\langle Q^2\rangle$")

ax.set_yscale("log")
ax.set_xscale("log")

ax.legend(loc="best")

ax.grid(alpha=.2)

fig.tight_layout()

fig.savefig(f"results/Q_var2.pdf", bbox_inches="tight")
plt.show()