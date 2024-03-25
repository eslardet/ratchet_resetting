import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from density import *
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)


alpha = 4.0
# beta_arr = [1.0]
beta = 1.0
gamma_arr = [0.2]
gamma = 0.2
t = format(10.0, '.6f')
# t_arr = np.arange(20, 21, 1)
# t_arr = [format(i, '.6f') for i in np.arange(0, 1.01, 0.1)] + [format(20, '.6f')]
t_arr = [format(i, '.6f') for i in [0.1, 0.5, 1.0, 5.0]]

# marker_list = ["o", "^", "s", "*", "D"]
marker_list = ["^", "s", "*", "D"]

legend_dict = {'d': r'$P_D(x)$', 'r': r'$P_R(x)$', 't': r'$P(x)$'}


phase = 'r'

save_plot = True

bins = 50
markersize=50


small = 22
big = 28

plt.rc('font', size=big)          # controls default text sizes
plt.rc('axes', labelsize=big)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=small)    # fontsize of the tick labels
plt.rc('ytick', labelsize=small)    # fontsize of the tick labels
plt.rc('legend', fontsize=small)    # legend fontsize

matplotlib.rcParams["font.family"] = 'STIXGeneral'
plt.rcParams['text.usetex'] = True

fig, ax = plt.subplots(figsize=(10, 6))

if phase == "d":
    ax.set_prop_cycle(color=plt.cm.Blues(np.linspace(0.2, 0.8, len(t_arr))))
elif phase == "r":
    ax.set_prop_cycle(color=plt.cm.Reds(np.linspace(0.2, 0.8, len(t_arr))))

for i, t in enumerate(t_arr):
    pos_file = get_pos_file(alpha, beta, gamma, t, phase)
    x_all = read_pos(pos_file)

    prob = get_prob_phase(alpha, beta, gamma, t, phase)
    # ax.hist(x_all, bins=100, density=True)
    hist, bin_edges = np.histogram(x_all, bins=bins, range=[0,alpha], density=True)
    bin_centers = 0.5*(bin_edges[1:]+bin_edges[:-1])

    # hist = np.append(hist, [hist[0]])
    # bin_centers = bin_edges[:-1]
    ax.scatter(bin_centers/alpha, prob*hist*alpha, marker=marker_list[i], s=markersize, label=r'$t={t}$'.format(t=str(float(t))))

    ax.plot(bin_centers/alpha, prob*hist*alpha, alpha=0.2, linewidth=3)

filename = 'P_{}_alpha{}_beta{}_gamma{}.csv'.format(phase, alpha, beta, gamma)
x_plot, y_plot = read_csv(filename)
ax.plot(x_plot/alpha, y_plot, color='k', linestyle="--", label = legend_dict[phase], zorder=6, linewidth=3)


ax.set_xlim(0,1)
ax.set_xlabel(r'$x$')
ax.set_ylabel(r'$P(x,t)$')
ax.legend(frameon=False)

if phase == "d":
    ax.set_ylim(bottom=0)
    ax.yaxis.set_minor_locator(MultipleLocator(0.1))
    ax.yaxis.set_major_locator(MultipleLocator(0.5))
else:
    ax.set_ylim(0,1)
    ax.yaxis.set_minor_locator(MultipleLocator(0.05))
    ax.yaxis.set_major_locator(MultipleLocator(0.2))

ax.xaxis.set_minor_locator(MultipleLocator(0.05))
ax.set_xticks([0, 0.25, 0.5, 0.75, 1], [r"$0$", r"$L/4$", r"$L/2$", r"$3L/4$", r"$L$"])


if save_plot:
    folder = os.path.abspath('./plots/density_time/')
    if not os.path.exists(folder):
        os.makedirs(folder)
    file_name = os.path.join(folder, phase + '_a' + str(alpha) + '_b' + str(beta) + '_c' + str(gamma) + '.pdf')
    plt.savefig(file_name, bbox_inches='tight')
    plt.close()
else:
    plt.show()