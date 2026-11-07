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

data = np.array([
    [100, 0.0200, 0.3963, 0.3964, 0.4356],
    [120, 0.0167, 0.4040, 0.4040, 0.1500],
    [140, 0.0143, 0.4101, 0.4099, 0.0506],
    [160, 0.0125, 0.4147, 0.4143, 0.0172],
    [180, 0.0111, 0.4189, 0.4183, 0.0057],
    [200, 0.0100, 0.4213, 0.4216, 0.0019],
    [220, 0.0091, 0.4253, 0.4253, 0.0006],
    [240, 0.0083, 0.4288, 0.4277, 0.0002],
    [260, 0.0077, 0.4306, 0.4315, 0.0001],
    [280, 0.0071, 0.4357, 0.4316, np.nan],
    [300, 0.0067, 0.4230, 0.4405, np.nan]
])

N_t = data[:, 0]
eta = data[:, 1]
d_star = data[:, 2]
d_star_dopo = data[:, 3]
attr = data[:, 4]

fig, ax = plt.subplots(figsize=(12,8))
ax.plot(data[:,0], data[:,2],'-o')

ax.set_xlabel(fr"$N_t$")
ax.set_ylabel(fr"$d^*$")

fig.tight_layout()

plt.show()
