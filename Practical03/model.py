# import aux
# import model
import warnings
# import numpy as np
# import matplotlib.pyplot as plt
# from scipy.integrate import odeint


def dGonorrhea(y, t, N, beta, partRate, rec):
    S, I = y
    # ODEs
    dSdt = -partRate * beta * I * S/N + rec*I
    dIdt = partRate * beta * I * S/N - rec*I
    return (dSdt, dIdt)


def dGonorrheaHeterogeneous(
            y, t, N,
            beta, cL, cMean, nH, nL, rec, gH, gL, cH
        ):
    SH, IH, SL, IL = y
    p = (gH*IH)/(nH*N)+(gL*IL)/(nL*N)
    # ODE
    dSHdt = -cH * beta * p * SH + rec * IH
    dIHdt = cH * beta * p * SH - rec * IH
    dSLdt = -cL * beta * p * SL + rec * IL
    dILdt = cL * beta * p * SL - rec * IL
    return (dSHdt, dIHdt, dSLdt, dILdt)


def dHIV(y, t, N, c, beta, gamma, mu, m, alpha):
    # Init
    S, I, A = y
    # ODE
    lmda = c * beta * I / (I + S)
    dSdt = ((alpha + m) * N) - (lmda * S + m * S)
    dIdt = (lmda * S) - (gamma * I + m * I)
    dAdt = (gamma * I) - ((m + mu) * A)
    return (dSdt, dIdt, dAdt)


def dHIVHeterogeneous(
            y, t, N,
            cH, cL, newH, newL, alpha,
            m, beta, gam, mu
        ):
    (SH, SL, IH, IL, AH, AL, DH, DL) = y
    # Intermediate params
    nH = IH/(IH+SH)
    nL = IL/(IL+SL)
    gH = (cH * nH) / (cH * nH + cL * nL)
    gL = 1 - gH
    # Partnership probability
    p = (gH * IH) / (SH + IH) + (gL * IL) / (SL + IL)
    # Rate changes
    (lambdaH, lambdaL) = (cH * beta * p,  cL * beta * p)
    dSH = newH * (alpha + m) * N - lambdaH * SH - m*SH
    dIH = lambdaH * SH - (gam + m) * IH
    dAH = gam * IH - (m + mu) * AH
    dDH = mu * AH
    dSL = newL * (alpha + m) * N - lambdaL * SL - m*SL
    dIL = lambdaL * SL - (gam + m) * IL
    dAL = gam * IL - (m + mu) * AL
    dDL = mu * AL
    return (dSH, dSL, dIH, dIL, dAH, dAL, dDH, dDL)
