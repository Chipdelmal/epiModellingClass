import aux
import model
import warnings
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
warnings.filterwarnings("ignore",category=DeprecationWarning)
plt.rcParams['figure.figsize'] = [12, 5]


def dGonorrhea(y, t, N, beta, partRate, rec):
    S, I = y
    # ODEs
    dSdt = -partRate*beta*I*S/N + rec*I
    dIdt = partRate*beta*I*S/N-rec*I
    return (dSdt, dIdt)


# User inputs
(N, tmax) = (200000000, 10)
(beta, partRate, duration) = (.75, 1, 2)
(iInit) = (1)
# Internals
(sInit, rec) = (N - iInit, 1 / duration)
t = np.linspace(0, tmax, tmax)
y0 = (sInit, iInit)
ret = odeint(dGonorrhea, y0, t, args=(N, beta, partRate, rec))
(S, I) = ret.T

I
2/12
