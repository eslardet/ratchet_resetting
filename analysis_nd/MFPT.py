import os
import csv
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)


small = 28
big = 30

plt.rc('font', size=40)          # controls default text sizes
plt.rc('axes', labelsize=big)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=small)    # fontsize of the tick labels
plt.rc('ytick', labelsize=small)    # fontsize of the tick labels
plt.rc('legend', fontsize=big)    # legend fontsize

matplotlib.rcParams["font.family"] = 'STIXGeneral'
plt.rcParams['text.usetex'] = True


# colors = plt.cm.viridis(np.linspace(0, 1, 3))
# colors = plt.cm.Reds(np.linspace(0.2, 1, len(d_list)))

l = 4.0
k = 4.0

d_list = [0.2, 0.5]
marker_list = ["^", "s"]

x_sim_plot = np.arange(0.0, 4.01, 0.2)
x_sim = [format(i, '.6f') for i in x_sim_plot]

fig, ax = plt.subplots(figsize=(10, 6))

def mfpt_analytic(l, k, d, x):
    if x < d*l:
        t = l*d/(k**2 * (np.exp(k) - 1)) * (l*(1-d)*np.exp(-k)*(np.exp(k*x/(l*d))-1) + k*x*(np.exp(k)-1) - l*(np.exp(k*x/(l*d))-1)*(1-k+d*(2*k-1)))
    else:
        t = l*(1-d)/(k**2 * (np.exp(k)-1)) * (l*d*(np.exp(-k*(x/l - d)/(1-d))-np.exp(-k)) + np.exp(k)*k*(l-x) - l*np.exp(k-k/(1-d)*(x/l - d)) * (d+k-2*d*k) + (x*k + l*d*(1-2*k)))

    return t

def get_sim_file(l, k, d, x):
    return './simulation_data_mfpt/ell{}_kappa{}_delta{}/mfpt_{}'.format(l, k, d, x)

for i, d in enumerate(d_list):
    x_plot = np.linspace(0.0, 4.0001, 100)
    t_analytic = []
    for x in x_plot:
        t_analytic.append(mfpt_analytic(l, k, d, x))
    ax.plot(x_plot/l, t_analytic)
    # ax.plot(k, eta, label=r'$\ell={}$'.format(l), color=colors[i])

    t_sim = []

    for x in x_sim:
        file_name = get_sim_file(l, k, d, x)

        with open(file_name) as f:
            reader = csv.reader(f, delimiter="\n")
            r = list(reader)
        t_mean = np.mean([float(i[0]) for i in r])

        t_sim.append(t_mean)
    
    ax.scatter(x_sim_plot/l, t_sim, s=100, marker=marker_list[i])



ax.set_xlabel(r'$\tilde{x}$', labelpad=10)
ax.set_ylabel(r'$\tilde{T}_R(\tilde{x})$', labelpad=10)

custom_lines = [plt.Line2D([0], [0], marker='^', markersize=10, color="tab:blue"),
                 plt.Line2D([0], [0], marker='s', markersize=10, color="tab:orange")]
labels = [r"$\delta=0.2$", r"$\delta=0.5$"]
ax.legend(custom_lines, labels, frameon=False)

ax.set_xlim(0, 1)
ax.set_ylim(bottom=0)


ax.xaxis.set_minor_locator(MultipleLocator(0.05))
ax.set_xticks([0, 0.25, 0.5, 0.75, 1], [r"$0$", r"$\ell/4$", r"$\ell/2$", r"$3\ell/4$", r"$\ell$"])

ax.yaxis.set_minor_locator(MultipleLocator(0.05))
ax.yaxis.set_major_locator(MultipleLocator(0.2))

ax.tick_params(axis='x', which='major', pad=5)
ax.tick_params(axis='y', which='major', pad=5)
ax.tick_params(which = 'both', direction='out', top=False, right=False)
ax.tick_params(which = 'major', length = 6)
ax.tick_params(which = 'minor', length = 3)




folder = os.path.abspath('./plots/mfpt/')
if not os.path.exists(folder):
    os.makedirs(folder)
filename = os.path.join(folder, 'MFPT_vs_x.pdf')
plt.savefig(filename, bbox_inches='tight')

# plt.show()
