import aux
import model
import warnings
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
warnings.filterwarnings("ignore", category=DeprecationWarning)
plt.rcParams['figure.figsize'] = [12, 5]

# User inputs #################################################################
(N, tmax) = (10000, 1000)
(beta, c) = (.05, 3)
(gamma, mu, m, alpha) = (
        1/9, 1/1,
        1/35, 0.025
    )
(iInit) = (1)
# Internals ###################################################################
(sInit, aInit) = ((N - iInit), 0)
t = np.linspace(0, tmax, tmax * 10)
y0 = (sInit, iInit, aInit)
# Run #########################################################################
ret = odeint(
        model.dHIV, y0, t, args=(N, c, beta, gamma, mu, m, alpha)
    )
(S, I, A) = ret.T
# Plot SIA ####################################################################
# tp = ((S, '#02146b', 'S'), (I, '#b4e830', 'I'), (A, '#ffb428', 'A'))
# (fig, ax) = aux.plotEpiDynamicsPop(tp, t/12, tmax/12, N)
# ax.set_xticks(np.arange(0, tmax/12, 25))
# ax.set_yticks(np.arange(0, N, 1000))
# plt.grid(
#         b=True, which='major', lw=.2, alpha=.5,
#         color='#666666', linestyle='--'
#     )
# Plot SIA ####################################################################
tp = (
        ((I+A)/N, '#02146b', 'Prevalence'),
        ((c*beta*I/(S+I)*S)/N, '#b4e830', 'Incidence')
    )
(fig, ax) = aux.plotEpiDynamics(tp, t/12, tmax/12, 1, ymax=.3)
ax.set_xticks(np.arange(0, tmax/12, 25))
ax.set_yticks(np.arange(0, 1, .1))
ax.set_ylim(0, 1)
plt.grid(
        b=True, which='major', lw=.2, alpha=.5,
        color='#666666', linestyle='--'
    )
