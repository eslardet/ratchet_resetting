import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import time
import os
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)

import sys
sys.path.append('./src_nd_python')
from ratchet_reset import Particle

seed = 20
np.random.seed(seed)

# Simulation parameters
ell = 4.0
kappa = 1.0
delta = 0.5
dt = 0.00001
x0 = 0

# samples = 10**3
total_t = 100

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


exp = np.exp
cos = np.cos
sin = np.sin
sqrt = np.sqrt
sinh = np.sinh
tanh = np.tanh
cosh = np.cosh

def P_D_Full(l,k,d):
    
    G_int = 2
    A = -k*(1 - exp(l))*(k + l*(d - 1))*(-d*l + k)*(d*l + k)*(exp(k/(1 - d)) - exp(k*(d - 2)/(d - 1)))*(-d*l + k + l)*exp(d*l + k/(d - 1))/(k*l*(-(-k + l*(d - 1))*(k + l*(2*d - 1))*(d*l - k)*exp(2*d*l) - (-k + l*(2*d - 1)*(k - 1))*(-d*l + k)*(-d*l + k + l)*exp(2*d*l + k) + (k + l*(d - 1))*(k + l*(2*d - 1)*(k - 1))*(d*l + k)*exp(k + l) - (k + l*(d - 1))*(d*l + k)*(-2*d*l + k + l)*exp(l)) + 2*(-k*(-d*l**4*(d - 1)*(k*(1 - 2*d)**2 - 1) - 2*k**4 + k**2*l**2*(2*d*(d - 1) + k*(1 - 2*d)**2 + 2) + (-d*l**4*(d - 1) + 2*k**4 - 2*k**2*l**2*(d*(d - 1) + 1))*cosh(k) + (d*l**4*(1 - 2*d)**2*(d - 1) + 2*k**4 - 2*k**2*l**2*(3*d*(d - 1) + 1))*sinh(k))*sinh(l/2) + l*(2*d**2*l**4*(d - 1)**2 + k**4 - k**2*l**2*(d*(d - 1) + k*(1 - 2*d)**2 + 1) + k**2*(-k**2 + l**2*(3*d*(d - 1) + 1))*sinh(k) - (2*d**2*l**4*(d - 1)**2 + k**4 + k**2*l**2*(-d**2 + d - 1))*cosh(k))*cosh(l/2))*exp(l*(d + 1/2)))
    
    return A*G_int

def energy_rate_analytic(l, k, d):
    # return (l*k*csch(l/2)*np.sinh(l*d/2)*np.sinh(l*(1-d)/2))/(d*(1-d))
    return (k*np.sinh(d*l/2)*np.sinh(l*(1-d)/2)) / (d*l*(1-d)*np.sinh(l/2)) * P_D_Full(l,k,d)


time_arr = []
total_energy_arr = []

t0 = time.time()

particle = Particle(x0, dt, ell, kappa, delta, seed)
for i in range(num_steps):

    time_arr.append(particle.t)
    total_energy_arr.append(particle.total_energy)

    particle.move()

fig, ax = plt.subplots(figsize=(10, 6))

ax.plot(time_arr, total_energy_arr, color='tab:red', label='Simulation')



print('Time taken: ', time.time() - t0)


time_plot = np.linspace(0, total_t, len(time_arr))
energy_analytic = energy_rate_analytic(ell, kappa, delta) * time_plot
ax.plot(time_plot, energy_analytic, color='black', linestyle='--', alpha=0.8, label='Analytic')


ax.set_xlabel(r'$\tilde{t}$')
ax.set_ylabel(r'$\tilde{E}_i$')
# ax.set_ylim(0, 1)
ax.set_xlim(0, total_t)
ax.legend(frameon=False)

# Tick locations
ax.xaxis.set_major_locator(MultipleLocator(5))
ax.xaxis.set_minor_locator(MultipleLocator(1))
ax.yaxis.set_major_locator(MultipleLocator(1))
ax.yaxis.set_minor_locator(MultipleLocator(0.25))


folder = 'plots/energy_rate/'
if not os.path.exists(folder):
    os.makedirs(folder)
# filename = 'L{}_a{}_h{}_D{}_r{}_seed{}'.format(L, a, h, D, r, seed)
filename = 'ell{}_kappa{}_delta{}_seed{}'.format(ell, kappa, delta, seed)
plt.savefig(folder + filename + '.pdf', bbox_inches='tight')

# plt.show()