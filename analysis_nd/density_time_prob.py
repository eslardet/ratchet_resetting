import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from density import *
import time

alpha = 4.0
alpha_arr = [1.0, 2.0, 3.0, 4.0]
# beta_arr = [1.0, 2.0, 3.0, 4.0]
beta = 1.0
# gamma_arr = np.round(np.arange(0.1, 0.51, 0.1),1)
gamma = 0.2
t = format(10.0, '.6f')
# t_arr = np.arange(20, 21, 1)
# t_arr = [format(i, '.6f') for i in np.arange(0, 1.01, 0.1)] + [format(20, '.6f')]
# t_arr = [format(i, '.6f') for i in [0.0, 0.1, 0.2, 0.4, 1.0, 20.0]]
t_arr = [format(i, '.6f') for i in np.arange(0,5.0,0.2)]


t0 = time.time()

small = 14
big = 20

plt.rc('font', size=small)          # controls default text sizes
plt.rc('axes', labelsize=big)     # fontsize of the x and y labels
plt.rc('xtick', labelsize=big)    # fontsize of the tick labels
plt.rc('ytick', labelsize=big)    # fontsize of the tick labels
plt.rc('legend', fontsize=small)    # legend fontsize

# matplotlib.rcParams["font.family"] = 'STIXGeneral'
# plt.rcParams['text.usetex'] = True


folder = os.path.abspath('./../plots/density_time/')
if not os.path.exists(folder):
    os.makedirs(folder)
file_name = 'beta' + str(beta) + "_gamma" + str(gamma)
save_file = open(os.path.join(folder, file_name + '.txt'), 'w')


fig, ax = plt.subplots(figsize=(7,5))


for alpha in alpha_arr:
    prob_arr = []
    save_file.write(alpha)
    for t in t_arr:
        prob = get_prob_phase(alpha, beta, gamma, t, phase='d')
        prob_arr.append(prob)

        save_file.write(t + "\t" + prob)

    t_arr_plot = [float(i) for i in t_arr]
    # ax.plot(t_arr_plot, prob_arr, '-o', label=r"$\alpha={}, \beta={}, \gamma={}$".format(alpha, beta, gamma))
    ax.plot(t_arr_plot, prob_arr, '-o', label=r"$\alpha={}$".format(alpha))

print("Time taken = " + str(time.time() - t0))

# ax.set_xlim(0,1)
ax.set_xlabel(r'$t$')
ax.set_ylabel(r'Probability in phase d$')
ax.legend()

# ax.set_ylim(0,1)

folder = os.path.abspath('./../plots/density_time/')
if not os.path.exists(folder):
    os.makedirs(folder)
file_name = os.path.join(folder, 'alpha' + str(alpha) + '.pdf')
plt.savefig(file_name, bbox_inches='tight')
plt.close()

save_file.close()