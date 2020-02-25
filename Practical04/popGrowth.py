
import model
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
plt.rcParams['figure.figsize'] = [12, 5]

# #############################################################################
# Init parameters
# #############################################################################
(N0, birthRate, simTime, reps) = (1, .5, 100, 100)

# #############################################################################
# Run the stochastic iterations of the Growth model
# #############################################################################
# t = np.linspace(0, simTime, simTime * 10)
# y0 = (N0, )
# ret = odeint(model.popGrowth, y0, t, args=tuple(birthRate))
traces = model.simulatePopGrowth(N0, birthRate, reps, simTime)

# #############################################################################
# Plot
# #############################################################################
fig = plt.figure(facecolor='w')
ax = fig.add_subplot(111, axisbelow=True)
for row in traces:
    ax.plot(range(simTime), row, "Red", alpha=0.5, lw=.1)
