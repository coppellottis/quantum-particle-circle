import numpy as np
import  matplotlib.pyplot as plt

def E(n, x) :
    return (n+x/2)**2


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

# x=theta/pi
x=np.linspace(-2,+2,1000)
y1= E(0,x)
y2 = E(1,x)
y3 = E(-1,x)

ax.plot(x,y1,'k:',alpha=.4,label=rf"$n=0$")
ax.plot(x,y2, 'k-.',alpha=.4,label=rf"$n=1$")
ax.plot(x,y3, 'k--',alpha=.4,label=rf"$n=-1$")

x4 = np.linspace(-2,-1,300)
y4 = E(1,x4)
x5 = np.linspace(-1,1,300) 
y5 = E(0,x5)
x6 = np.linspace(1,2,300)
y6 = E(-1,x6)
ax.plot(x4,y4,'r',lw=2)
ax.plot(x5,y5,'r',lw=2)
ax.plot(x6,y6,'r',lw=2,label="Stato fondamentale")

ax.set_xlabel(rf'$\theta/\pi$')
ax.set_ylabel(rf"$E_n(\theta)/(\hbar^2/2mR^2)$")

ax.legend(loc="best")

ax.grid(alpha=.2)
ax.set_xlim(-2,2)

fig.tight_layout()

fig.savefig(f"results/energy_levels.pdf", bbox_inches="tight")
plt.show()