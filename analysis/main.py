import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from density import *

h = 1
a = 0.2
# a_arr = np.round(np.arange(0.1,0.6,0.1),1)
a_arr = [0.2]
D = 0.1
r = 0.2
r_arr = [0.9]
# t = 20
t = format(1.0, '.6f')
t_arr = np.arange(20, 21, 1)
# t_arr = [format(i, '.6f') for i in np.arange(0, 1.01, 0.1)] + [20]
# t_arr = [format(1.0, '.6f')]
phase = 'r'

plot_pos_phase(h, a_arr, D, r_arr, t_arr, phase, save_plot=False, show_plot=True, bins=100)
# plot_pos_both(h, a, D, r, t, save_plot=False, show_plot=True)
# plot_pos_total(h, a, D, r, t_arr)