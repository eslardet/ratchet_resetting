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

legend_dict = {'d': r'$\tilde{P}_D(\tilde{x})$', 'r': r'$\tilde{P}_R(\tilde{x})$', 't': r'$\tilde{P}(\tilde{x})$'}

save_plot = True

bins = 40
markersize = 100


small = 28
big = 30

plt.rc('font', size=40)          # controls default text sizes
plt.rc('axes', labelsize=big)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=small)    # fontsize of the tick labels
plt.rc('ytick', labelsize=small)    # fontsize of the tick labels
plt.rc('legend', fontsize=big)    # legend fontsize

matplotlib.rcParams["font.family"] = 'STIXGeneral'
plt.rcParams['text.usetex'] = True

labels = ['(a)', '(b)', '(c)']
phases = ["d", "r"]

fig, axs = plt.subplots(1,3, figsize=(42,8))

for j in range(2):
    ax = axs[j]
    phase = phases[j]
    if phase == "d":
        ax.set_prop_cycle(color=plt.cm.Blues(np.linspace(0.2, 0.8, len(t_arr))))
    elif phase == "r":
        ax.set_prop_cycle(color=plt.cm.Reds(np.linspace(0.2, 0.8, len(t_arr))))

    filename = 'P_{}_alpha{}_beta{}_gamma{}.csv'.format(phase, alpha, beta, gamma)
    x_plot, y_plot = read_csv(filename)
    # ax.plot(x_plot/alpha, y_plot, color='k', linestyle="--", label = legend_dict[phase], linewidth=5, zorder=4)
    ax.plot(x_plot/alpha, y_plot, color='k', linestyle="--", label = legend_dict[phase], linewidth=5, zorder=4)

    for i, t in enumerate(t_arr):
        pos_file = get_pos_file(alpha, beta, gamma, t, phase)
        x_all = read_pos(pos_file)

        prob = get_prob_phase(alpha, beta, gamma, t, phase)
        # ax.hist(x_all, bins=100, density=True)
        hist, bin_edges = np.histogram(x_all, bins=bins, range=[0,alpha], density=True)
        bin_centers = 0.5*(bin_edges[1:]+bin_edges[:-1])

        # boundary_hist = np.histogram(x_all, bins=100, range=[0,alpha], density=True)[0][0]
        # hist = np.array([boundary_hist] + hist.tolist() + [boundary_hist])
        # bin_centers = np.array([0] + bin_centers.tolist() + [alpha])

        
        # bin_centers = np.append(bin_centers, np.array([alpha], dtype=float))
        # hist = np.append(hist, [boundary_hist])

        # bin_centers = np.append(np.array([0], dtype=float), bin_centers)
        # hist = np.append(np.array([boundary_hist]), hist)

        # hist = np.append(hist, [hist[0]])
        # bin_centers = bin_edges[:-1]

        # if t != "5.000000":
        # ax.plot(bin_centers/alpha, prob*hist*alpha, alpha=0.5, linewidth=5, zorder=2)
        ax.plot(bin_centers/alpha, prob*hist, alpha=0.5, linewidth=5, zorder=2)

        # ax.scatter(bin_centers/alpha, prob*hist*alpha, marker=marker_list[i], s=markersize, label=r'$\tilde{t}={t}$'.format(t=str(float(t))), zorder=3)
        ax.scatter(bin_centers/alpha, prob*hist, marker=marker_list[i], s=markersize, label=r'$\tilde{t}=' + str(float(t)) + r"$", zorder=3)

        



    ax.text(-0.12, 0.96, labels[j], horizontalalignment='center', verticalalignment='center', transform=ax.transAxes)

    ax.set_xlim(0,1)
    ax.set_xlabel(r'$\tilde{x}$')
    if phase == "d":
        ax.set_ylabel(r'$\tilde{P}_D(\tilde{x},\tilde{t})$', labelpad=20)
    elif phase == "r":
        ax.set_ylabel(r'$\tilde{P}_R(\tilde{x},\tilde{t})$', labelpad=20)
    ax.legend(frameon=False)

    if phase == "d":
        ax.set_ylim(bottom=0)
        ax.yaxis.set_minor_locator(MultipleLocator(0.05))
        ax.yaxis.set_major_locator(MultipleLocator(0.1))
    else:
        ax.set_ylim(0,0.3)
        # ax.set_ylim(bottom=0)
        ax.yaxis.set_minor_locator(MultipleLocator(0.05))
        ax.yaxis.set_major_locator(MultipleLocator(0.1))

    ax.xaxis.set_minor_locator(MultipleLocator(0.05))
    ax.set_xticks([0, 0.25, 0.5, 0.75, 1], [r"$0$", r"$\ell/4$", r"$\ell/2$", r"$3\ell/4$", r"$\ell$"])


ax = axs[2]
ax.text(-0.12, 0.96, labels[2], horizontalalignment='center', verticalalignment='center', transform=ax.transAxes)

folder = os.path.abspath('./plots/density_time_prob/')
if not os.path.exists(folder):
    os.makedirs(folder)
file_name = 'beta' + str(beta) + "_gamma" + str(gamma)
save_file = os.path.join(folder, file_name + '.txt')

with open(save_file) as f:
    reader = csv.reader(f, delimiter="\n")
    r = list(reader)

num_alpha = int(len(r)/3) 

ax.set_prop_cycle(color=plt.cm.BuPu(np.linspace(0.2, 1, num_alpha)))

for k in range(num_alpha):
    alpha = float(r[3*k][0])

    t_arr = r[3*k+1][0].split('\t')[:-1]
    t_plot = [float(i) for i in t_arr]
    prob = r[3*k+2][0].split('\t')[:-1]
    prob_plot = [float(i) for i in prob]
    ax.plot(t_plot, prob_plot, "-o", label=r"$\ell=" + str(int(alpha)) + r"$", linewidth=5, markersize=10)


ax.set_xlim(0,5)
ax.set_ylim(0.5, 1.0)
ax.set_xlabel(r'$\tilde{t}$')
ax.set_ylabel(r'$p_D(\tilde{t})$', labelpad=10)
ax.legend(frameon=True, loc="upper right")


ax.xaxis.set_minor_locator(MultipleLocator(0.25))
ax.yaxis.set_minor_locator(MultipleLocator(0.05))


if save_plot:
    folder = os.path.abspath('./plots/density_time/')
    if not os.path.exists(folder):
        os.makedirs(folder)
    # file_name = os.path.join(folder, 'both_a' + str(alpha) + '_b' + str(beta) + '_c' + str(gamma) + '.pdf')
    file_name = os.path.join(folder, "DensityVsTime_ProbabilityDiffusion.pdf")
    plt.savefig(file_name, bbox_inches='tight')
    plt.close()
else:
    plt.show()