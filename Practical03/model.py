

def dGonorrhea(y, t, N, beta, partRate, rec):
    S, I = y
    # ODEs
    dSdt = -partRate*beta*I*S/N + rec*I
    dIdt = partRate*beta*I*S/N-rec*I
    return (dSdt, dIdt)


def dGonorrheaHeterogeneous(y, t, N, beta, cl, cm, nh, nl, rec):
    Sl, Sh, Il, Ih = y
    # Auxiliary functions
    ch = (cm - cl * nl) / (nh)
    gh = (ch * nh) / (ch * nh + cl * nl)
    gl = 1 - gh
    pt = gh * Ih + gl * Il
    ll = cl * beta * pt
    lh = ch * beta * pt
    # ODEs
    dSldt = -ll * Sl + rec * Il
    dShdt = -lh * Sh + rec * Ih
    dIldt = ll * Sl - rec * Il
    dIhdt = lh * Sh - rec * Ih
    return (dSldt, dShdt, dIldt, dIhdt)
