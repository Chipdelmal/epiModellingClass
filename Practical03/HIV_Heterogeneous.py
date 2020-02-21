import aux
import model
import warnings
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
warnings.filterwarnings("ignore", category=DeprecationWarning)
plt.rcParams['figure.figsize'] = [12, 5]

(N, Inf0, tmax) = (10000, 1, 100)

(cH, cL) = (8, 0.2)
(alpha, beta, newH) = (0.025, .05, 0.15)
(gam, mu, m) = (1/9, 1/1, 1/35)
newL = 1 - newH
NH = newH * N

(SH, SL, IH, IL, AH, AL, DH, DL) = (
        NH - Inf0, N - NH,
        Inf0, 0,
        0, 0,
        0, 0
    )

t = np.linspace(0, tmax, tmax * 2)
y0 = (SH, SL, IH, IL, AH, AL, DH, DL)
ret = odeint(
        model.dHIVHeterogeneous, y0, t,
        args=(
                N, cH, cL, newH, newL, alpha,
                m, beta, gam, mu
            )
    )
(SH, SL, IH, IL, AH, AL, DH, DL) = ret.T
# Plot SIA ####################################################################
Nt = SH + IH + AH + SL + IL + AL
PREV = (IH + IL + AH + AL) / Nt
PREV_H = (IH + AH) / (SH + IH + AH)
PREV_L = (IL + AL) / (SL + IL + AL)
# Plot SIA ####################################################################
tp = (
        (PREV,      '#02146b', 'Prevalence'),
        (PREV_H,    '#b4e830', 'Prevalence high'),
        (PREV_L,    '#e21e7b', 'Prevalence IL')
    )
(fig, ax) = aux.plotEpiDynamics(tp, t/12, tmax/12, 1, ymax=.3)
ax.set_xticks(np.arange(0, tmax/12, 25))
ax.set_yticks(np.arange(0, 1, .1))
ax.set_ylim(0, 1)
plt.grid(b=True, which='major', lw=.2, alpha=.5,
         color='#666666', linestyle='--')

(PREV[-1], PREV_H[-1], PREV_L[-1])
