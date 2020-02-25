
import model
import numpy as np
import matplotlib.pyplot as plt
# from scipy.integrate import odeint
plt.rcParams['figure.figsize'] = [12, 5]

# #############################################################################
# Init parameters
# #############################################################################
(N0, birthRate, deathRate, simTime, reps) = (20, .5, .3, 200, 1000)

# #############################################################################
# Run the stochastic iterations of the Growth model
# #############################################################################
traces = model.simulatePopGrowthDeath(N0, birthRate, deathRate, reps, simTime)

# #############################################################################
# Plot
# #############################################################################
fig = plt.figure(facecolor='w')
ax = fig.add_subplot(111, axisbelow=True)
for row in traces:
    ax.plot(range(simTime), row, "Red", alpha=0.5, lw=.05)
