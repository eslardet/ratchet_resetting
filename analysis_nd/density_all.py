import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from density import *


alpha = 4.0
beta = 1.0
gamma = 0.5
t = format(20.0, '.6f')
bins = 100

fig, ax = initalize_plotting()

 
for phase in ['d', 'r']:
    filename = 'P_{}_alpha{}_beta{}_gamma{}.csv'.format(phase, alpha, beta, gamma)
    x_plot, y_plot = read_csv(filename)
    ax.plot(x_plot/alpha, y_plot, label="Phase " + phase)

filename = 'P_t_alpha{}_beta{}_gamma{}.csv'.format(alpha, beta, gamma)
x_plot, y_plot = read_csv(filename)
ax.plot(x_plot/alpha, y_plot, label="Total", color='black')

# x_plot, y_plot = read_csv(filename)
# ax.plot(x_plot, y_plot, label='Total')

x_all = []
    
for phase in ['d', 'r']:
    pos_file = get_pos_file(alpha, beta, gamma, t, phase)
    x = read_pos(pos_file)
    x_all += list(x)
    
    prob = get_prob_phase(alpha, beta, gamma, t, phase)
    # ax.hist(x_all, bins=100, density=True)
    hist, bin_edges = np.histogram(x, bins=bins, range=[-0.01, 4.01], density=True)
    bin_centers = 0.5*(bin_edges[1:]+bin_edges[:-1])/alpha
    ax.scatter(bin_centers, hist*prob*alpha, label="Phase " + phase)

hist, bin_edges = np.histogram(x_all, bins=bins, density=True)
bin_centers = 0.5*(bin_edges[1:]+bin_edges[:-1])/alpha
ax.scatter(bin_centers, hist*alpha, label='Total', color='black')



# ax.set_xlim(0,1)
ax.set_xticks([0, 0.25, 0.5, 0.75, 1], [r"$0$", r"$L/4$", r"$L/2$", r"$3L/4$", r"$L$"])
ax.set_xlabel(r'$x$')
ax.set_ylabel(r'$P(x)$')
ax.legend()

plt.show()