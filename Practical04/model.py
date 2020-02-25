
import math
import numpy as np


# #############################################################################
# Simple population growth model
# #############################################################################
def popGrowth(y, t, N0, birthRate):
    pop = y
    dPopdt = birthRate * pop
    return (dPopdt, )


def simPopGrowth(simTime, N0, birthRate):
    # Init Array
    pop = np.zeros(simTime)
    pop[0] = int(N0)
    # Iterate through stochastic traces
    for t in range(1, simTime):
        currPop = int(pop[t-1])
        births = np.random.binomial(1, birthRate, size=currPop)[0]
        pop[t] = pop[t-1] + births
    return pop


def simulatePopGrowth(N0, birthRate, reps, simTime):
    traces = np.empty((reps, simTime))
    for rep in range(reps):
        traces[rep] = simPopGrowth(simTime, N0, birthRate)
    return traces


# #############################################################################
# Population Growth Model with Death
# #############################################################################
def simPopGrowthDeath(simTime, N0, birthRate, deathRate):
    pop = np.zeros(simTime)
    pop[0] = int(N0)
    # Iterate through stochastic traces
    for t in range(1, simTime):
        currPop = int(pop[t-1])
        actions = np.random.binomial(1, birthRate + deathRate, size=currPop)
        bd = np.random.multinomial(
                1, [birthRate, deathRate],
                size=np.sum(actions)
            ).T
        (births, deaths) = [np.sum(i) for i in bd]
        pop[t] = pop[t-1] + births - deaths
    return pop


def simulatePopGrowthDeath(N0, birthRate, deathRate, reps, simTime):
    traces = np.empty((reps, simTime))
    for rep in range(reps):
        traces[rep] = simPopGrowthDeath(simTime, N0, birthRate, deathRate)
    return traces


# #############################################################################
# SIR model
# #############################################################################
def sir(u, parms, t):
    bet, gamm, iota, N, dt = parms
    S, I, R, Y = u
    lambd = bet*(I+iota)/N
    ifrac = 1.0 - math.exp(-lambd*dt)
    rfrac = 1.0 - math.exp(-gamm*dt)
    infection = np.random.binomial(S, ifrac)
    recovery = np.random.binomial(I, rfrac)
    return [S-infection, I+infection-recovery, R+recovery, Y+infection]


def simulateSIR(params, tf, tl, iZero):
    t = np.linspace(0, tf, tl)
    (S, I, R, Y) = [np.zeros(tl) for i in range(4)]
    u = [params[3] - iZero, iZero, 0, 0]
    S[0], I[0], R[0], Y[0] = u
    for j in range(1, tl):
        u = sir(u, params, t[j])
        S[j], I[j], R[j], Y[j] = u
    return {'t': t, 'S': S, 'I': I, 'R': R, 'Y': Y}
