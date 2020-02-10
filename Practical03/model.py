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
