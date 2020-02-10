import aux
# import model
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


def dGonorrheaHeterogeneous(y, t, N, beta, cL, cMean, nH, nL, rec):
    SH, IH, SL, IL = y
    p = (gH*IH)/(nH*N)+(gL*IL)/(nL*N)
    # ODE
    dSHdt = -cH*beta*p*SH+rec*IH
    dIHdt = cH*beta*p*SH-rec*IH
    dSLdt = -cL*beta*p*SL+rec*IL
    dILdt = cL*beta*p*SL-rec*IL
    return (dSHdt, dIHdt, dSLdt, dILdt)


############################################################
(N, tmax)  = (200000000, 15 * 12)
(beta, rec) = (3/4, 1/2)
(nH, cMean, cL)  = (0.02, 2/12, 1.4/12)
Inf0 = 1
############################################################
nL = 1-nH
cH = (cMean-cL*nL)/nH
gH = cH*nH/(cH*nH+cL*nL)
gL = 1-gH
############################################################
SH0 = N*nH-Inf0
IH0 = Inf0
SL0 = N-N*nH
IL0 = 0
y0 = (SH0, IH0, SL0, IL0)
t = np.linspace(0, tmax, tmax * 10)
ret = odeint(
        dGonorrheaHeterogeneous, y0, t,
        args=(N, beta, cL, cMean, nH, nL, rec)
    )
(SH, IH, SL, IL) = ret.T
############################################################
prev = (IH+IL)/N
prevH = IH/(N*nH)
prevL = IL/(N*nL)
tp = (
        (SH, '#02146b', 'SH'), (IH, '#ffb428', 'IH'),
        (SL, '#b4e830', 'SL'), (IL, '#e21e7b', 'IL')
    )
(fig, ax) = aux.plotEpiDynamicsPop(tp, t/12, tmax/12, 1)
# ax.set_xticks(np.arange(0, tmax/12, 1))
# ax.set_yticks(np.arange(0, 1, .1))
ax.set_ylim(0, np.amax(SH))
ax.set_xlim(0, tmax / 12)
plt.grid(
        b=True, which='major', lw=.2, alpha=.5,
        color='#666666', linestyle='--'
    )

tp = (
        (prev, '#ffb428', 'Prevalence'),
        (prevH, '#02146b', 'Prevalence IH'),
        (prevL, '#e21e7b', 'Prevalence IL')
    )
(fig, ax) = aux.plotEpiDynamicsPop(tp, t/12, tmax/12, 1)
# ax.set_xticks(np.arange(0, tmax/12, 1))
# ax.set_yticks(np.arange(0, 1, .1))
ax.set_ylim(0, .3)
ax.set_xlim(0, tmax / 12)
plt.grid(
        b=True, which='major', lw=.2, alpha=.5,
        color='#666666', linestyle='--'
    )
