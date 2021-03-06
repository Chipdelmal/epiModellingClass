
import matplotlib.pyplot as plt


def plotEpiDynamics(triplets, t, tmax, N, alpha=.5, lw=4, ymax=1):
    fig = plt.figure(facecolor='w')
    ax = fig.add_subplot(111, axisbelow=True)
    for tp in triplets:
        ax.plot(t, tp[0]/N, tp[1], alpha=0.5, lw=4, label=tp[2])
    ax.set_xlabel('Time (days)')
    ax.set_ylabel('Population Fraction')
    ax.set_ylim(0, ymax)
    ax.set_xlim(0, tmax)
    ax.legend()
    return fig, ax


def plotEpiDynamicsPop(triplets, t, tmax, N, alpha=.5, lw=4):
    fig = plt.figure(facecolor='w')
    ax = fig.add_subplot(111, axisbelow=True)
    for tp in triplets:
        ax.plot(t, tp[0], tp[1], alpha=0.5, lw=lw, label=tp[2])
    ax.set_xlabel('Time (years)')
    ax.set_ylabel('Population Fraction')
    ax.set_ylim(0, N)
    ax.set_xlim(0, tmax)
    ax.legend()
    return fig, ax
