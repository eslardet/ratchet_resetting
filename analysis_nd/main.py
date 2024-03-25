import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from density import *

alpha = 4.0
beta_arr = [1.0]
gamma_arr = [0.2]
t = format(10.0, '.6f')
# t_arr = np.arange(20, 21, 1)
t_arr = [format(i, '.6f') for i in np.arange(0, 1.01, 0.1)] + [format(20, '.6f')]

# t_arr = [format(10.0, '.6f')] + [20]
# t_arr = [format(10.0, '.6f')]
phase = 'd'

# plot_pos_phase(h, a_arr, D, r_arr, t_arr, phase, save_plot=False, show_plot=True, bins=50)
# plot_pos_both(h, a, D, r, t, save_plot=False, show_plot=True)
# plot_pos_total(h, a, D, r, t_arr)


plot_pos_phase(alpha, beta_arr, gamma_arr, t_arr, phase, save_plot=False, show_plot=True, bins=50)

# fig, ax = initalize_plotting()
# # ax = plot_compare(h, a, D, r, t, phase, bins=50, ax=ax)
# filename = 'P_{}_h{}_a{}_D{}_r{}.csv'.format(phase, h, a, D, r)
# ax = plot_csv(filename, ax=ax)
# ax.set_xlim(0,1)
# ax.set_xlabel(r'$x$')
# ax.set_ylabel(r'$P(x)$')

# folder = os.path.abspath('./plots/density/')
# if not os.path.exists(folder):
#     os.makedirs(folder)
# filename = os.path.join(folder, 'compare_h' + str(h) + '_a' + str(a) + '_D' + str(D) + '_r' + str(r) + '_t' + str(t) + '.pdf')
# # plt.savefig(filename, bbox_inches='tight')

# plt.show()

