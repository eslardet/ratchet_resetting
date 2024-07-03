import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import time
import os
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)

import sys
sys.path.append('./src_nd_python')
from ratchet_reset import Particle

seed = 1
np.random.seed(seed)

# Simulation parameters
alpha = 1.0
beta = 0.2
gamma = 0.05
dt = 0.001
x0 = 0

# samples = 10**3
total_t = 2.5

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
plt.rcParams['axes.labelpad']=10

phase_colors = {'d': 'tab:blue', 'r': 'tab:red'}

fig, ax = plt.subplots(figsize=(10, 6))



final_loc = []

x_loc = []
last_t = 0
previous_x = x0
current_phase = 'd'
t0 = time.time()
# for s in range(samples):
# Create particle
particle = Particle(x0, dt, alpha, beta, gamma, seed)
for i in range(num_steps):

    dx = np.abs(particle.x - previous_x)
    if particle.phase == current_phase and dx < 0.5*alpha:
        x_loc.append(particle.x)
        previous_x = particle.x
    else:
        current_phase = particle.phase
        # t_plot = np.arange(last_t, particle.t, dt)
        t_plot = np.linspace(last_t, particle.t, len(x_loc))
        x_loc = np.array(x_loc)/alpha
        ax.plot(t_plot, x_loc, color=phase_colors[current_phase])
        x_loc = [particle.x]
        last_t = particle.t
        previous_x = particle.x

    particle.move()

final_loc.append(particle.x)


t_plot = np.linspace(last_t, particle.t, len(x_loc))
x_loc = np.array(x_loc)/alpha
ax.plot(t_plot, x_loc, color=phase_colors[current_phase])

print('Time taken: ', time.time() - t0)


# ax.hist(final_loc, bins=100)
# t_plot = np.arange(0, total_t, dt)
# ax.plot(t_plot, x_loc)



ax.set_xlabel(r'$\tilde{t}$')
ax.set_ylabel(r'$\tilde{x}/\alpha$')
ax.set_ylim(0, 1)
ax.set_xlim(0, total_t)

# Tick locations


folder = 'plots/trajectory/'
if not os.path.exists(folder):
    os.makedirs(folder)
# filename = 'L{}_a{}_h{}_D{}_r{}_seed{}'.format(L, a, h, D, r, seed)
filename = 'alpha{}_beta{}_gamma{}_seed{}_unwrapped'.format(alpha, beta, gamma, seed)
# plt.savefig(folder + filename + '.pdf', bbox_inches='tight')

plt.show()