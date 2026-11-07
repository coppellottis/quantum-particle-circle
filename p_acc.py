from scipy.special import erfc, erfcinv
import numpy as np
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


# Nt = 100
data_100 = np.array([
    [1.0251, 0.9061],
    [1.0507, 0.8079],
    [1.0771, 0.7152],
    [1.1049, 0.6266],
    [1.1349, 0.5371],
    [1.1683, 0.4525],
    [1.2073, 0.3632],
    [1.2563, 0.2676],
    [1.4340, 0.0861],
    [1.3920, 0.1132],
    [1.3290, 0.1713]
])

# Nt = 50
data_50 = np.array([
    [1.0355, 0.9201],
    [1.0717, 0.8440],
    [1.1090, 0.7679],
    [1.1483, 0.6981],
    [1.1908, 0.6258],
    [1.2380, 0.5555],
    [1.2931, 0.4836],
    [1.3625, 0.4069],
    [1.4652, 0.3114],
    [1.5544, 0.2516],
    [1.6138, 0.2226]
])

data_200 = np.array([
    [1.0178, 0.9001],
    [1.0358, 0.8015],
    [1.0545, 0.7092],
    [1.0742, 0.6115],
    [1.0954, 0.5223],
    [1.1190, 0.4273],
    [1.1466, 0.3340],
    [1.1812, 0.2398],
    [1.2326, 0.1395],
    [1.2772, 0.0839],
    [1.3069, 0.0580]
])

x = np.linspace(0,6)
y = erfc(x/(2*np.sqrt(2)))

fig, ax = plt.subplots(figsize=(12,8))

ax.plot(x,y,'k-.', label=r"erfc$((C-1)\sqrt{N_t}/2\sqrt{2})$")

ax.scatter((data_50[:,0]-1)*np.sqrt(50),data_50[:,1], label=r"$N_t=50$",s=50, marker='o')
ax.scatter((data_100[:,0]-1)*np.sqrt(100),data_100[:,1], label=r"$N_t=100$",s=50, marker='o')
ax.scatter((data_200[:,0]-1)*np.sqrt(200),data_200[:,1], label=r"$N_t=200$",s=50, marker='o')

ax.set_xlabel(r"$(C-1)\sqrt{N_t}$")
ax.set_ylabel("Accettazione dello scambio")

ax.grid(alpha=.2)
ax.legend(loc="best")

fig.tight_layout()

fig.savefig("results/p_acc.pdf", bbox_inches="tight")

plt.show()