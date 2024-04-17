import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from density import *
import time
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)

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

small = 22
big = 28

plt.rc('font', size=big)          # controls default text sizes
plt.rc('axes', labelsize=big)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=small)    # fontsize of the tick labels
plt.rc('ytick', labelsize=small)    # fontsize of the tick labels
plt.rc('legend', fontsize=small)    # legend fontsize

matplotlib.rcParams["font.family"] = 'STIXGeneral'
plt.rcParams['text.usetex'] = True




folder = os.path.abspath('./plots/density_time_prob/')
if not os.path.exists(folder):
    os.makedirs(folder)
file_name = 'beta' + str(beta) + "_gamma" + str(gamma)
save_file = os.path.join(folder, file_name + '.txt')


fig, ax = plt.subplots(figsize=(10, 6))


with open(save_file) as f:
    reader = csv.reader(f, delimiter="\n")
    r = list(reader)

num_alpha = int(len(r)/3) 

ax.set_prop_cycle(color=plt.cm.Blues(np.linspace(0.2, 1, num_alpha)))

for k in range(num_alpha):
    alpha = float(r[3*k][0])

    t_arr = r[3*k+1][0].split('\t')[:-1]
    t_plot = [float(i) for i in t_arr]
    prob = r[3*k+2][0].split('\t')[:-1]
    prob_plot = [float(i) for i in prob]
    ax.plot(t_plot, prob_plot, "-o", label=r"$\ell=" + str(int(alpha)) + r"$")


ax.set_xlim(0,5)
ax.set_ylim(0.5, 1.0)
ax.set_xlabel(r'$t$')
ax.set_ylabel(r'$p_D(t)$', labelpad=10)
ax.legend(frameon=True, loc="upper right")


ax.xaxis.set_minor_locator(MultipleLocator(0.25))
ax.yaxis.set_minor_locator(MultipleLocator(0.05))

# ax.set_ylim(0,1)


plt.show()

# plt.savefig(os.path.join(folder, file_name + '.pdf'), bbox_inches='tight')
# plt.close()
