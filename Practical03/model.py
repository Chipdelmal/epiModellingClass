

def dGonorrhea(y, t, N, beta, partRate, rec):
    S, I = y
    # ODEs
    dSdt = -partRate*beta*I*S/N + rec*I
    dIdt = partRate*beta*I*S/N-rec*I
    return (dSdt, dIdt)
