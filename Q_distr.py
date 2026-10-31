import numpy as np
import  matplotlib.pyplot as plt

def P(Q, beta) :
    Q1 = np.arange(-50,+51)
    den = np.sum(np.exp(-Q1**2/(2*beta)))
    return np.exp(-Q**2/(2*beta))/den


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

beta = 2
x = np.arange(-10,11)

betas = [0.5, 1, 2]

for j, beta in enumerate(betas):
    arr = []
    width = .2

    for i in x:
        y = P(i, beta)
        arr.append([i, y])

    arr = np.array(arr)

    np.set_printoptions(suppress=True, precision=3)
    print(f"beta={beta}: \n{np.column_stack([arr[:,0], arr[:,1]*100])}\n")

    ax.bar(arr[:,0]+width*(1-j),arr[:,1],width=width, label=fr"$\beta={beta}$",alpha=.6)

ax.set_xlabel(rf'$Q$')
ax.set_ylabel(rf"$P(Q)$")

ax.set_xlim(-5,5)
ax.set_xticks(np.arange(-5, 6, 1))

ax.legend(loc="best")

ax.grid(alpha=.2)

fig.tight_layout()

fig.savefig(f"results/Q_distr.pdf", bbox_inches="tight")
plt.show()