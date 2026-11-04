import numpy as np
from scipy.special import erf
import  matplotlib.pyplot as plt
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

def P(Q, eta, N_t):
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


eps = 1e-10
theta = 0
beta = 2

Nt = [10, 12, 50, 100]

for n in Nt:
    print(f"Nt: {n}, <Q^2>: {Q_var(beta/n,n)}, <Q^4>: {Q_2var(beta/n,n)}, b2: {b2(beta/n,n)}\n")

#fig, ax = plt.subplots(figsize=(12, 8))

# Nt = np.arange(10,102,2)
# yy = []

# for n in Nt:
#     yy.append([n, b2(2/n,n)])

# y = np.asarray(yy)

# ax.plot(y[:,0],y[:,1],"-o")

# ax.set_yscale("log")
# ax.set_xlabel(fr"$N_t$")
# ax.set_ylabel(fr"$b_2$")

# ax.grid(alpha=.2)

# for j, N in enumerate(Nt):

#     n = np.arange(0, 21)
#     y = []

#     for i in n:
#         y.append([i, eigenvalue(i,theta,beta/N)])

#     y = np.asarray(y)

#     ax.plot(y[:,0],abs(y[:,1]/y[0,1]),"-o", label=fr"$N_t={N},\;\eta={beta/N}$")
#     mask = y[:,1]<0
#     ax.plot(y[mask,0], abs(y[mask,1]/y[0,1]), "*", color="k", linestyle="None", markersize=10)

# ax.set_yscale("log")
# ax.set_xlabel(fr"$n$")
# ax.set_ylabel(r"$|\lambda_n(0)/\lambda_0(0)|$")

# ax.set_xticks(np.arange(0, 21, 2))

# ax.legend(loc="best")
# ax.grid(alpha=.2)

# fig.savefig(f"results/discrete.pdf", bbox_inches="tight")

# fig2, ax2 = plt.subplots(figsize=(12,8))

# NN_t = np.arange(10,101,2)
# yy = []
# for k in NN_t:
#     yy.append([k, Q_var(2/k, k), Q_2var(2/k,k)])

# yy = np.asarray(yy)
# ax2.plot(yy[:,0],yy[:,1],'-o',label=fr"$\langle Q^2\rangle$")
# ax2.plot(yy[:,0],yy[:,2],'-o',label=fr"$\langle Q^4\rangle$")
# ax2.plot(yy[:,0],3*(yy[:,1]**2),'-o',label=fr"$3\langle Q^2\rangle^2$")
# ax2.axhline(y=2,color='black', linestyle="-.",label=r"$\langle Q^2\rangle=\hat\beta=2$")

# ax2.set_yscale("log")
# ax2.grid(alpha=.2)

# ax2.set_xlabel(fr"$N_t$")
# #ax2.set_ylabel(fr"$\langle Q\rangle$")
# ax2.legend(loc="best")

# fig2.savefig(f"results/discrete4.pdf",bbox_inches="tight")

# plt.show()