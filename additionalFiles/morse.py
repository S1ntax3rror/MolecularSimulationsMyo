
import numpy as np
import matplotlib.pyplot as plt

r1 = np.arange(1.0, 4.0, 0.1)
r2 = np.arange(1.03, 4.03, 0.1)

rre = ????
rDe = -????
rbeta = ???

def morse(d1, d2, re, De, beta):
    V1 = -De*(1. - np.exp(-beta*(d1 - re)))**2
    V2 = -De*(1. - np.exp(-beta*(d2 -   re)))**2
    return V1 + V2

V = []
for d1, d2 in zip(r1, r2):
    V.append(morse(d1, d2, rre, rDe, rbeta))

V = np.array(V)
rng = np.random.default_rng()
V_noise = 0.02 * rng.normal(size=V.size)
V += V_noise

#print(r1.shape, r2.shape, )
l = np.concatenate(
        (r1.reshape(1, -1), r2.reshape(1, -1), np.array(V).reshape(1, -1)),
        axis=0)

plt.plot(r1, V)
plt.show()


np.savetxt("dimorse_test.dat", l.T)
