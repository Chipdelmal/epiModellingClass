def dEbola(y, t, N, betaI, betaD, f, m, r, mu):
    S, E, I, R, D, B = y
    # ODEs
    dSdt = -(betaI * I / N + betaD * D/N) * S
    dEdt = (betaI * I / N + betaD * D/N) * S - f * E
    dIdt = f * E - (m + r) * I
    dRdt = r * I
    dDdt = m * I - mu * D
    dBdt = mu * D
    return dSdt, dEdt, dIdt, dRdt, dDdt, dBdt


def dEbolaBurial(y, t, N, betaI, betaD, f, m, r, mu, actTime=90):
    S, E, I, R, D, B = y
    # Interventions Conditions
    if t > actTime:
        mu = mu * 2
    # ODEs
    dSdt = -(betaI * I / N + betaD * D/N) * S
    dEdt = (betaI * I / N + betaD * D/N) * S - f * E
    dIdt = f * E - (m + r) * I
    dRdt = r * I
    dDdt = m * I - mu * D
    dBdt = mu * D
    return dSdt, dEdt, dIdt, dRdt, dDdt, dBdt


def dEbolaCRBodies(y, t, N, betaI, betaD, f, m, r, mu, actTime=90):
    S, E, I, R, D, B = y
    # Interventions Conditions
    if t > actTime:
        betaD = betaD / 2
    # ODEs
    dSdt = -(betaI * I / N + betaD * D/N) * S
    dEdt = (betaI * I / N + betaD * D/N) * S - f * E
    dIdt = f * E - (m + r) * I
    dRdt = r * I
    dDdt = m * I - mu * D
    dBdt = mu * D
    return dSdt, dEdt, dIdt, dRdt, dDdt, dBdt


def dEbolaCRInfected(y, t, N, betaI, betaD, f, m, r, mu, actTime=90):
    S, E, I, R, D, B = y
    # Interventions Conditions
    if t > actTime:
        betaI = betaI / 2
    # ODEs
    dSdt = -(betaI * I / N + betaD * D/N) * S
    dEdt = (betaI * I / N + betaD * D/N) * S - f * E
    dIdt = f * E - (m + r) * I
    dRdt = r * I
    dDdt = m * I - mu * D
    dBdt = mu * D
    return dSdt, dEdt, dIdt, dRdt, dDdt, dBdt


def dEbolaCRAll(y, t, N, betaI, betaD, f, m, r, mu, actTime=90):
    S, E, I, R, D, B = y
    # Interventions Conditions
    if t > actTime:
        betaI = betaI / 2
        betaD = betaD / 2
        mu = mu * 2
    # ODEs
    dSdt = -(betaI * I / N + betaD * D/N) * S
    dEdt = (betaI * I / N + betaD * D/N) * S - f * E
    dIdt = f * E - (m + r) * I
    dRdt = r * I
    dDdt = m * I - mu * D
    dBdt = mu * D
    return dSdt, dEdt, dIdt, dRdt, dDdt, dBdt
