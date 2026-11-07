from scipy.special import erfc, erfcinv
import numpy as np

Pacc = .6
Nt = 400
C = 1+2*np.sqrt(2)*erfcinv(Pacc)/np.sqrt(Nt)

print(f"C: {C}")
print(f"erfcinv(.4): {erfcinv(Pacc)}")

beta0 = 2
eta0 = beta0/Nt

for i in np.arange(0,10):
    eta = (C**i)*eta0
    beta = eta*Nt
    print(f"Replica: {i}, eta: {eta}, beta: {beta}")


xx = [.9,.8,.7,.6,.5,.4,.3,.2,.1,.05,.03]
for x in xx:
    C = 1+2*np.sqrt(2)*erfcinv(x)/np.sqrt(Nt)
    print(f"Pacc={x}, C={C}")

beta_max = 4
eta_max = beta_max/Nt

NR = np.sqrt(Nt)/(2*np.sqrt(2)*erfcinv(Pacc))*np.log(eta_max/eta0)+1
print(f"NR per beta_max = {beta_max}: {NR}")