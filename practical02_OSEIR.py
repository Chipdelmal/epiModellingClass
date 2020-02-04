import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt


###############################################################################
# Functions
###############################################################################
# The SEIR model differential equations.
def deriv(y, t, N, betaI, betaH, f, m, r, mu):
    S, E, I, R, D, B = y
    dSdt = -(betaI * I / N + betaH * D/N) * S
    dEdt = (betaI * I / N + betaH * D/N) * S - f * E
    dIdt = f * E - (m + r) * I
    dRdt = r * I
    dDdt = m * I - mu * D
    dBdt = mu * D
    return dSdt, dEdt, dIdt, dRdt, dDdt, dBdt


###############################################################################
# Main
###############################################################################
(N, tmax) = (1000, 365)
(betaI, betaH) = (0.3, 1)
(f, r, m, mu) = (1/6, 1/10, 1/7.5, 1)
(E0, I0, R0, D0, B0) = (0, 1, 0, 0, 0)
S0 = N - (E0 + I0 + R0 + D0 + B0)

t = np.linspace(0, tmax, tmax)
y0 = S0, E0, I0, R0, D0, B0
ret = odeint(deriv, y0, t, args=(N, betaI, betaH, f, m, r, mu))
(S, E, I, R, D, B) = ret.T

###############################################################################
# Plots
###############################################################################
# Plot the data on three separate curves for S(t), I(t) and R(t)
fig = plt.figure(facecolor='w')
ax = fig.add_subplot(111, axisbelow=True)
ax.plot(t, S/N, 'b', alpha=0.5, lw=2, label='Susceptible')
ax.plot(t, E/N, 'm', alpha=0.5, lw=2, label='Exposed')
ax.plot(t, I/N, 'r', alpha=0.5, lw=2, label='Infected')
ax.plot(t, R/N, 'g', alpha=0.5, lw=2, label='Recovered')
ax.plot(t, B/N, 'r', alpha=0.5, lw=2, label='Buried')
ax.plot(t, D/N, 'g', alpha=0.5, lw=2, label='Dead')
ax.set_xlabel('Time /days')
ax.set_ylabel('Number (1000s)')
ax.set_ylim(0, 1)
ax.yaxis.set_tick_params(length=0)
ax.xaxis.set_tick_params(length=0)
ax.grid(b=True, which='major', c='w', lw=2, ls='-')
legend = ax.legend()
legend.get_frame().set_alpha(0.5)
plt.show()
