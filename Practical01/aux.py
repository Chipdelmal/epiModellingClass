import matplotlib.pyplot as plt


def plotEpiDynamics(triplets, t, N, alpha=.5):
    fig = plt.figure(facecolor='w')
    ax = fig.add_subplot(111, axisbelow=True)
    for tp in triplets:
        ax.plot(t, tp[0]/N, tp[1], alpha=0.5, lw=2, label=tp[2])
    ax.set_xlabel('Time (days)')
    ax.set_ylabel('Population Fraction')
    ax.set_ylim(0, 1)
    ax.set_xlim(0, len(t))
    ax.legend()
    return fig, ax