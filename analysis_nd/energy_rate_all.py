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

pi = np.pi
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

def get_sim_file(l, k, d, t, s=None):
    
    if s is not None:
        filepath = './simulation_data/ell{}_kappa{}_delta{}/s{}/energy_{}'.format(l, k, d, s, t)
    else:
        filepath = './simulation_data/ell{}_kappa{}_delta{}/energy_{}'.format(l, k, d, t)
    return filepath

def energy_rate_sim(l, k, d, t, seed=None):
    file_name = get_sim_file(l, k, d, t, seed)
    with open(file_name) as f:
        reader = csv.reader(f, delimiter="\n")
        r = list(reader)
    energy_all = np.array([float(i[0]) for i in r])

    # file_name = get_sim_file(l, k, d, t_ss)
    # with open(file_name) as f:
    #     reader = csv.reader(f, delimiter="\n")
    #     r = list(reader)
    # energy_ss = np.array([float(i[0]) for i in r])

    # energy_all = energy_all - energy_ss

    # t_total = float(t) - float(t_ss)
    t_total = float(t)

    energy_mean = np.mean(energy_all / t_total)
    energy_std = np.std(energy_all / t_total)
    return energy_mean, energy_std


l = 4.0
k_arr = [1.0, 2.0, 3.0, 4.0]
# k_arr = [2.0, 3.0]
d_arr = np.linspace(0.01, 0.99, 100)
# d_arr_sim = np.round(np.arange(0.1, 0.91, 0.1),1)
d_arr_sim = [0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95]

t = format(11000000.0, '.6f')
# t_ss = format(20.0, '.6f')

# seed = None
seed = 2

colors = ["tab:blue", "tab:orange", "tab:green", "tab:red"]
marker_list = ["^", "s", "*", "D"]

fig, axs = plt.subplots(2,2, figsize=(20,16))

# for i, k in enumerate(k_arr):
for i, ax in enumerate(fig.axes):
    # ax = axs[i]
    k = k_arr[i]
    all_E = []
    all_E_std = []
    for d in d_arr:  
        all_E.append(energy_rate_analytic(l, k, d))
    ax.plot(d_arr, all_E, label=r'$\kappa = {}$'.format(k), color=colors[i])


    all_E = []
    all_E_std = []
    for d in d_arr_sim:
        E_mean, E_std = energy_rate_sim(l, k, d, t, seed)
        all_E.append(E_mean)
        all_E_std.append(E_std)
    all_E = np.array(all_E)

    print(all_E)
    ax.scatter(d_arr_sim, all_E, color=colors[i], marker=marker_list[i], s=200)

    # ax.errorbar(d_arr_sim, all_E, yerr=all_E_std, fmt=marker_list[i], color=colors[i], label=r'$\kappa = {}$'.format(k), capsize=5, markersize=10)

    # ax.legend(frameon=False)
    ax.set_xlabel(r'$\delta$', labelpad=10)
    ax.set_ylabel(r'$\dot{E}_i$', labelpad=10)

    ax.set_xlim(0, 1)
    ax.xaxis.set_minor_locator(MultipleLocator(0.05))
    ax.xaxis.set_major_locator(MultipleLocator(0.2))

# ax = axs[0,0]
# ax.yaxis.set_minor_locator(MultipleLocator(0.005))
# ax.yaxis.set_major_locator(MultipleLocator(0.01))

# ax = axs[1,1]
# ax.yaxis.set_minor_locator(MultipleLocator(0.005))
# ax.yaxis.set_major_locator(MultipleLocator(0.02))


custom_lines = [plt.Line2D([0], [0], marker=marker_list[i], markersize=10, color=colors[i]) for i in range(len(k_arr))]
labels = [r"$\kappa=" + str(k) + "$" for k in k_arr]
fig.legend(custom_lines, labels, frameon=False, loc="upper center", ncol=4, bbox_to_anchor=(0.5, 1.1))

plt.tight_layout()


folder = os.path.abspath('./plots/energy_rate/')
if not os.path.exists(folder):
    os.makedirs(folder)
file_name = 'energy_rate_l4_'
if seed is not None:
    file_name += 'seed' + str(seed) + '_'
file_name += 'dt10-5_t' + str(round(float(t)))
save_file = os.path.join(folder, file_name + '.png')
plt.savefig(save_file, bbox_inches='tight')

# plt.show()

