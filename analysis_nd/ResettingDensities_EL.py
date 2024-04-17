import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
from density import *


small = 28
big = 30

plt.rc('font', size=40)          # controls default text sizes
plt.rc('axes', labelsize=big)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=small)    # fontsize of the tick labels
plt.rc('ytick', labelsize=small)    # fontsize of the tick labels
plt.rc('legend', fontsize=big)    # legend fontsize
markersize = 100

# matplotlib.rcParams["font.family"] = "serif"
plt.rcParams['text.usetex'] = True
# matplotlib.rcParams['mathtext.fontset'] = 'cm'
matplotlib.rcParams['font.family'] = 'STIXGeneral'



labels = ['(a)', '(b)', '(c)', '(d)']
alpha_list = [4.0,4.0,4.0,1.0]
beta_list = [1.0,1.0,1.0,4.0]
gamma_list_analytical = ["0.00001",0.2,0.5,"0.00001"]
gamma_list = [0.0,0.2,0.5,0.0]
bins = 20
t = format(20.0, '.6f')

color_dict = {'d': 'tab:blue', 'r': 'tab:red', 't': 'black'}
symbol_dict = {'d': 'o', 'r': 's', 't': '^'}
legend_dict = {'d': r'$P_D(x)$', 'r': r'$P_R(x)$', 't': r'$P(x)$'}


fig, axs = plt.subplots(2,2, figsize=(24,16))


for i, ax in enumerate(fig.axes):
    alpha = alpha_list[i]
    beta = beta_list[i]
    gamma = gamma_list[i]
    gamma_analytical = gamma_list_analytical[i]

    ## Analytical ##
    for phase in ['t', 'd', 'r']:
        filename = 'P_{}_alpha{}_beta{}_gamma{}.csv'.format(phase, alpha, beta, gamma_analytical)
        x_plot, y_plot = read_csv(filename)
        ax.plot(x_plot/alpha, y_plot, color=color_dict[phase], label = legend_dict[phase])

    ## Numerical ##
    x_all = []

    for phase in ['d', 'r']:
        pos_file = get_pos_file(alpha, beta, gamma, t, phase)
        x = read_pos(pos_file)
        x_all += list(x)

    hist, bin_edges = np.histogram(x_all, bins=bins, density=True)
    bin_centers = 0.5*(bin_edges[1:]+bin_edges[:-1])/alpha
    ax.scatter(bin_centers, hist*alpha, color=color_dict['t'], s=markersize, marker=symbol_dict['t'], label = legend_dict['t'] + " (simulation)")

    for phase in ['d', 'r']:
        pos_file = get_pos_file(alpha, beta, gamma, t, phase)
        x = read_pos(pos_file)
        
        prob = get_prob_phase(alpha, beta, gamma, t, phase)
        # ax.hist(x_all, bins=100, density=True)
        hist, bin_edges = np.histogram(x, bins=bins, density=True)
        bin_centers = 0.5*(bin_edges[1:]+bin_edges[:-1])/alpha
        ax.scatter(bin_centers, hist*prob*alpha, color=color_dict[phase], s=markersize, marker=symbol_dict[phase], label = legend_dict[phase] + " (simulation)")


    ## General plotting ##
    ax.text(-0.12, 0.96, labels[i], horizontalalignment='center', verticalalignment='center', transform=ax.transAxes)
    
    ax.xaxis.set_minor_locator(MultipleLocator(0.05))
    ax.set_xticks([0, 0.25, 0.5, 0.75, 1], [r"$0$", r"$L/4$", r"$L/2$", r"$3L/4$", r"$L$"])
    ax.yaxis.set_minor_locator(MultipleLocator(0.05))
    ax.yaxis.set_major_locator(MultipleLocator(0.2))

    ax.tick_params(axis='x', which='major', pad=5)
    ax.tick_params(axis='y', which='major', pad=5)
    ax.tick_params(which = 'both', direction='out', top=False, right=False)
    ax.tick_params(which = 'major', length = 6)
    ax.tick_params(which = 'minor', length = 3)
    
    ax.set_xlabel(r'$x$')
    ax.set_ylabel(r'$P(x), P_D(x), P_R(x)$', labelpad=10)

    ax.set_xlim(0, 1)
    ax.set_ylim(bottom=0)

    ax.vlines(gamma, 0, 1.4, linestyle = 'dashed', alpha = 0.1, linewidth = 3.5, color = 'black')

    if i == 3:
        ax.legend(frameon=False, loc='center', ncol=2, handletextpad=0.2, columnspacing = 0.7)

plt.tight_layout()

folder = os.path.abspath('./plots/density/')
if not os.path.exists(folder):
    os.makedirs(folder)
plt.savefig(folder + '/density_alpha{}_beta{}_gamma{}.pdf'.format(alpha, beta, gamma), bbox_inches='tight')

# plt.show()
