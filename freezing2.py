import numpy as np
import emcee
import pandas as pd
import  matplotlib.pyplot as plt
from scipy.special import erf
import math
import cmath

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

def eigenvalue(n, theta, eta):
    k = theta + 2*math.pi*n

    return (np.sqrt(2*math.pi*eta)
            * np.exp(-eta*k**2/2)
            * np.real(erf(1/(2*np.sqrt(2*eta))
                          - 1j*k*np.sqrt(eta/2))))

def Z(theta, eta, Nt):
    Z = eigenvalue(0, theta, eta)**Nt
    for n in range(1, 10):
        pair = (eigenvalue(n, theta, eta)**Nt
                + eigenvalue(-n, theta, eta)**Nt)

        Z += pair

        if abs(pair) / abs(Z) < eps:
            break
    return Z

def P2(Q, eta, N_t):
    N_theta = 1000
    theta = np.linspace(-math.pi,math.pi, N_theta)
    P_Q = 0
    for th in theta:
        P_Q += np.exp(1j*th*Q)*(Z(th,eta,N_t)/Z(0,eta,N_t))
    return P_Q/N_theta

def Q_var(eta,N_t):
    Q = np.arange(-10,11)
    var = 0
    for q in Q:
        var+= (q**2)*P(q,eta,N_t)
    return var

def Q_2var(eta,N_t):
    Q = np.arange(-10,11)
    var = 0
    for q in Q:
        var+= (q**4)*P(q,eta,N_t)
    return var

def Q_mean(eta,N_t):
    Q = np.arange(-10,11)
    mean = 0
    for q in Q:
        mean += q*P(q,eta,N_t)
    return mean

def b2(eta,N_t):
    return (Q_2var(eta,N_t)-3*(Q_var(eta,N_t)**2))/(12*Q_var(eta,N_t))

def P(Q, beta) :
    Q1 = np.arange(-50,+51)
    den = np.sum(np.exp(-Q1**2/(2*beta)))
    return np.exp(-Q**2/(2*beta))/den

N = 10
n = N
beta = 2
eps = 1e-10
theta = 0

fig, ax = plt.subplots(figsize=(12, 8))
width = .2

eta = beta/N
data = pd.read_csv(f"data/top_freezing2/N{N}_a{eta:.4f}_beta{beta:.4f}.csv")

QQ = np.arange(-6,7)
hist = [np.sum(data["Q"] == i)/np.size(data,axis=0) for i in QQ]
ax.bar(QQ+width,hist,width=width,label=fr"$N_t={n}$, MC (locale)")

yQ = np.array([P(Q, beta) for Q in QQ])
    
ax.bar(QQ, yQ,width=width,alpha=.4,label=fr"Teoria continua")

yQ2 = np.array([P2(Q, beta/n,n) for Q in QQ])

ax.bar(QQ-width,yQ2,width=width,alpha=.4,label=fr"Teoria discreta")
ax.legend(loc="best")

fig.tight_layout()

fig.savefig("results/freezing2b.pdf",bbox_inches="tight")


# N = np.arange(1,5)*100
# beta = 2
# eps = 1e-10
# theta = 0

# fig, ax = plt.subplots(4,1,figsize=(8,12),sharex=True)

# width = .2

# for i, n in enumerate(N):
#     eta = beta/n
#     data = pd.read_csv(f"data/top_freezing/N{n}_a{eta:.4f}_beta{beta:.4f}.csv")

#     QQ = np.arange(-6,7)
#     hist = [np.sum(data["Q"] == i)/np.size(data,axis=0) for i in QQ]
#     ax.flat[i%4].bar(QQ+width,hist,width=width,label=fr"$N_t={n}$")

#     yQ = np.array([P(Q, beta) for Q in QQ])
    
#     ax.flat[i%4].bar(QQ, yQ,width=width,alpha=.4)

#     yQ2 = np.array([P2(Q, beta/n,n) for Q in QQ])

#     ax.flat[i%4].bar(QQ-width,yQ2,width=width,alpha=.4)
#     ax.flat[i%4].legend(loc="best")

# ax[-1].set_xlabel("Q")
# fig.tight_layout()

# fig.savefig("results/freezing2.pdf",bbox_inches="tight")


plt.show()