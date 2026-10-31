import numpy as np
import  matplotlib.pyplot as plt
import math

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

def chi1(theta, beta):
    return 0.5*beta*(theta[:,None]**2)

def chi2(theta,beta):
    return 2*np.exp(-1/(2*beta))*(1-np.cos(theta[:,None]))

x = np.linspace(-1,1,1000)

fig, ax = plt.subplots(figsize=(12, 8))

y1 = chi1(x*math.pi,2)/chi1(np.array([math.pi]),2)
ax.plot(x,y1,label=r"$\beta=2$")

y2 = chi2(x*math.pi,0.3)/chi2(np.array([math.pi]),0.3)
ax.plot(x,y2,label=r"$\beta=0.3$")

ax.set_xlabel(rf'$\theta/\pi$')
ax.set_ylabel(rf"$(f(\theta)-f(0))/(f(\pi)-f(0))$")

ax.legend(loc="best")

ax.grid(alpha=.2)

fig.tight_layout()

fig.savefig(f"results/free_energy.pdf", bbox_inches="tight")
plt.show()