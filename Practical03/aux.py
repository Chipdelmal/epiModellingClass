
import numpy as np
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


def parseEbolaPlotTriplets(S, E, I, R, B, D):
    tp = (
            (S, '#02146b', 'S'), (E, '#ffb428', 'E'), (I, '#b4e830', 'I'),
            (R, '#e21e7b', 'R'), (B, '#888888', 'B'), (D, '#12eaea', 'D')
        )
    return tp
