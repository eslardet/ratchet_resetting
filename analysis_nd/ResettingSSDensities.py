import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
from density import *


small = 24
big = 30

plt.rc('font', size=big)          # controls default text sizes
plt.rc('axes', labelsize=big)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=small)    # fontsize of the tick labels
plt.rc('ytick', labelsize=small)    # fontsize of the tick labels
plt.rc('legend', fontsize=small)    # legend fontsize
markersize = 100

# matplotlib.rcParams["font.family"] = "serif"
plt.rcParams['text.usetex'] = True
# matplotlib.rcParams['mathtext.fontset'] = 'cm'
matplotlib.rcParams['font.family'] = 'STIXGeneral'



labels = ['(a)', '(b)', '(c)', '(d)']
# ell_list = [4.0,4.0,4.0,1.0]
# kappa_list = [1.0,1.0,1.0,4.0]
# delta_list_analytical = ["0.00001",0.2,0.5,"0.00001"]
# delta_list = [0.0,0.2,0.5,0.0]

ell = 4.0
kappa = 1.0
delta = 0.2

bins = 20
t = format(20.0, '.6f')

color_dict = {'d': 'blue', 'r': 'red', 't': 'black'}
symbol_dict = {'d': 'o', 'r': 's', 't': '^'}
legend_dict = {'d': r'$\tilde{P}_D(\tilde{x})$', 'r': r'$\tilde{P}_R(\tilde{x})$', 't': r'$\tilde{P}(\tilde{x})$'}


# fig, axs = plt.subplots(2,2, figsize=(24,16))
fig, ax = plt.subplots(figsize=(10,6))


# for i, ax in enumerate(fig.axes):
#     ell = ell_list[i]
#     kappa = kappa_list[i]
#     delta = delta_list[i]
#     delta_analytical = delta_list_analytical[i]

## Analytical ##
for phase in ['t', 'd', 'r']:
    filename = 'P_{}_ell{}_kappa{}_delta{}.csv'.format(phase, ell, kappa, delta)
    x_plot, y_plot = read_csv(filename, model="v2")
    ax.plot(x_plot/ell, y_plot, color=color_dict[phase], label = legend_dict[phase])

## Numerical ##
x_all = []

for phase in ['d', 'r']:
    pos_file = get_pos_file(ell, kappa, delta, t, phase, model="v2")
    x = read_pos(pos_file)
    x_all += list(x)

hist, bin_edges = np.histogram(x_all, bins=bins, density=True)
bin_centers = 0.5*(bin_edges[1:]+bin_edges[:-1])/ell
ax.scatter(bin_centers, hist, color=color_dict['t'], s=markersize, marker=symbol_dict['t'], label = legend_dict['t'] + " (simulation)")

for phase in ['d', 'r']:
    pos_file = get_pos_file(ell, kappa, delta, t, phase, model="v2")
    x = read_pos(pos_file)
    
    prob = get_prob_phase(ell, kappa, delta, t, phase, model="v2")
    # ax.hist(x_all, bins=100, density=True)
    hist, bin_edges = np.histogram(x, bins=bins, density=True)
    bin_centers = 0.5*(bin_edges[1:]+bin_edges[:-1])/ell
    ax.scatter(bin_centers, hist*prob, color=color_dict[phase], s=markersize, marker=symbol_dict[phase], label = legend_dict[phase] + " (simulation)")


## General plotting ##
# ax.text(-0.12, 0.96, labels[i], horizontalalignment='center', verticalalignment='center', transform=ax.transAxes)

ax.xaxis.set_minor_locator(MultipleLocator(0.05))
ax.set_xticks([0, 0.25, 0.5, 0.75, 1], [r"$0$", r"$\ell/4$", r"$\ell/2$", r"$3\ell/4$", r"$\ell$"])
ax.yaxis.set_minor_locator(MultipleLocator(0.05))
ax.yaxis.set_major_locator(MultipleLocator(0.2))

ax.tick_params(axis='x', which='major', pad=5)
ax.tick_params(axis='y', which='major', pad=5)
ax.tick_params(which = 'both', direction='out', top=False, right=False)
ax.tick_params(which = 'major', length = 6)
ax.tick_params(which = 'minor', length = 3)

ax.set_xlabel(r'$\tilde{x}$')
ax.set_ylabel(r'$\tilde{P}(\tilde{x}), \tilde{P}_D(\tilde{x}), \tilde{P}_R(\tilde{x})$', labelpad=10)

ax.set_xlim(0, 1)
ax.set_ylim(0, 1)


ax.vlines(delta, 0, 1.4, linestyle = 'dotted', alpha = 0.1, linewidth = 3.5, color = 'black')

# if i == 3:
ax.legend(frameon=False, loc='upper center', ncol=2, handletextpad=0.2, columnspacing = 0.7)

plt.tight_layout()

folder = os.path.abspath('./plots/density_stop_start/')
if not os.path.exists(folder):
    os.makedirs(folder)
plt.savefig(folder + '/density_ell{}_kappa{}_delta{}_sim.pdf'.format(ell, kappa, delta), bbox_inches='tight')

# plt.show()
