import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt


###############################################################################
# Functions
###############################################################################
# The SEIR model differential equations.
def deriv(y, t, N, beta, f, r):
    S, E, I, R = y
    dSdt = -beta * S * I / N
    dEdt = beta * S * I / N - f * E
    dIdt = f * E - r * I
    dRdt = r * I
    return dSdt, dEdt, dIdt, dRdt


###############################################################################
# Main
###############################################################################
# Total population, N.
N = 1000
latentPeriod = 8
infectiousPeriod = 7
R_0 = 3.25
R0 = 0
I0 = 1
E0 = 0
# Everyone else, S0, is susceptible to infection initially.
S0 = N - I0 - R0 - E0
# Contact rate, beta, and mean recovery rate, gamma, (in 1/days).
(f, r) = (1/latentPeriod, 1/infectiousPeriod)
beta = R_0 * r
# A grid of time points (in days)
t = np.linspace(0, 160, 160)


# Initial conditions vector
y0 = S0, E0, I0, R0
# Integrate the SIR equations over the time grid, t.
ret = odeint(deriv, y0, t, args=(N, beta, f, r))
S, E, I, R = ret.T

###############################################################################
# Plots
###############################################################################
# Plot the data on three separate curves for S(t), I(t) and R(t)
fig = plt.figure(facecolor='w')
ax = fig.add_subplot(111, axisbelow=True)
ax.plot(t, S/N, 'b', alpha=0.5, lw=2, label='Susceptible')
ax.plot(t, E/N, 'm', alpha=0.5, lw=2, label='Exposed')
ax.plot(t, I/N, 'r', alpha=0.5, lw=2, label='Infected')
ax.plot(t, R/N, 'g', alpha=0.5, lw=2, label='Recovered with immunity')
ax.set_xlabel('Time /days')
ax.set_ylabel('Number (1000s)')
ax.set_ylim(0, 1)
ax.yaxis.set_tick_params(length=0)
ax.xaxis.set_tick_params(length=0)
ax.grid(b=True, which='major', c='w', lw=2, ls='-')
legend = ax.legend()
legend.get_frame().set_alpha(0.5)
plt.show()
