import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import time
import os
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)

import sys
sys.path.append('./src')
from ratchet_reset import Particle

seed = 1
np.random.seed(seed)

# Simulation parameters
L = 1
a = 0.2
h = 0.1
D = 0.1
dt = 0.001
r = 2.0
x0 = 0

# samples = 10**3
total_t = 5.0

# num_steps = 20000
num_steps = int(total_t/dt)

small = 22
big = 28

plt.rc('font', size=big)          # controls default text sizes
plt.rc('axes', labelsize=big)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=small)    # fontsize of the tick labels
plt.rc('ytick', labelsize=small)    # fontsize of the tick labels
plt.rc('legend', fontsize=small)    # legend fontsize

matplotlib.rcParams["font.family"] = "serif"
plt.rcParams['text.usetex'] = True
plt.rcParams['axes.labelpad']=5

phase_colors = {'d': 'tab:blue', 'r': 'tab:red'}

fig, ax = plt.subplots(figsize=(8,10))



final_loc = []

x_loc = []
last_t = 0
previous_x = x0
current_phase = 'd'
t0 = time.time()
# for s in range(samples):
# Create particle
particle = Particle(x0, dt, D, L, a, h, r, seed)
for i in range(num_steps):

    dx = np.abs(particle.x - previous_x)
    if particle.phase == current_phase and dx < 0.5*L:
        x_loc.append(particle.x)
        previous_x = particle.x
    else:
        # t_plot = np.arange(last_t, particle.t, dt)
        if dx < 0.5*L:
            x_loc.append(particle.x)
            print(particle.x)
            print(particle.t)
        elif L-previous_x < 0.1:
            x_loc.append(L)
        elif previous_x < 0.1:
            x_loc.append(0)
        t_plot = np.linspace(last_t, particle.t, len(x_loc))
        ax.plot(x_loc, t_plot, color=phase_colors[current_phase])
        x_loc = [particle.x]
        last_t = particle.t
        previous_x = particle.x
        current_phase = particle.phase
    particle.move()

final_loc.append(particle.x)


t_plot = np.linspace(last_t, particle.t, len(x_loc))
ax.plot(x_loc, t_plot, color=phase_colors[current_phase])

print('Time taken: ', time.time() - t0)


# ax.hist(final_loc, bins=100)
# t_plot = np.arange(0, total_t, dt)
# ax.plot(t_plot, x_loc)

ax.vlines(a, 0, total_t, color='black', linestyle='--', alpha=0.2, linewidth=5)

ax.set_ylabel(r'$t$')
ax.set_xlabel(r'$x$')
ax.set_xlim(0, L)
ax.set_ylim(0, total_t)
# ax.set_ylim(0.265, total_t)


## Make custom ticks
# ax.set_yticks([0,5])
# ax.set_xticks(ax.get_xticks(), [0, '', '', '', '', ''])
ax.set_xticks([0,0.5,1])
ax.set_xticks(ax.get_xticks(), [0, r"$L/2$", r'$L$'])
# ax.set_xticks([0,0.2,0.4,0.6,0.8,1])
# ax.set_xticks(ax.get_xticks(), [0, r'$L/5$', r"$2L/5$", r"$3L/5$", r"$4L/5$", r'$L$'])

ax.xaxis.set_minor_locator(MultipleLocator(0.1))
ax.yaxis.set_major_locator(MultipleLocator(1))
ax.yaxis.set_minor_locator(MultipleLocator(0.5))

# plt.tick_params(left=False, right=False)

ax.plot((0), (1.01), ls="", marker="^", ms=10, color="k", transform=ax.get_xaxis_transform(), clip_on=False)

folder = os.path.abspath('./plots/trajectory/')
if not os.path.exists(folder):
    os.makedirs(folder)
filename = 'L{}_a{}_h{}_D{}_r{}_seed{}_inverse'.format(L, a, h, D, r, seed)
plt.savefig(os.path.join(folder, filename + '.pdf'), bbox_inches='tight')
plt.savefig(os.path.join(folder, filename + '.svg'), bbox_inches='tight')
# plt.savefig(os.path.join(folder, filename + '.png'), bbox_inches='tight')
# print(os.path.join(folder, filename + '.png'))
# plt.show()