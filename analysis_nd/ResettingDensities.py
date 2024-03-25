# -*- coding: utf-8 -*-
"""
Created on Fri Jan  7 14:09:41 2022

@author: crobe
"""

import time
import datetime
start_time = time.time()
print('Program started at', datetime.datetime.now().time())

import numpy as np
import csv
import matplotlib.pyplot as plt
import matplotlib
import pickle
from matplotlib.lines import Line2D
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
matplotlib.rcParams['mathtext.fontset'] = 'cm'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
to_rgba = matplotlib.colors.to_rgba
pi = np.pi
exp = np.exp
cos = np.cos
sin = np.sin
sqrt = np.sqrt

###############################################################################

############################# load in data files ##############################

###############################################################################

# subfigure a

x_list_a = []
density_a = []
diff_density_a = []
reset_density_a = []

# full density
with open("P_delta=0_ell=4_kappa=1.csv", "r") as file:
   lines = csv.reader(file)
   for line in lines:
       x = float(line[0])
       P = float(line[1])
       x_list_a.append(x)
       density_a.append(P)

# diffusion phase density    
with open("PD_delta=0_ell=4_kappa=1.csv", "r") as file:
   lines = csv.reader(file)
   for line in lines:
       x = float(line[0])
       P = float(line[1])
       diff_density_a.append(P)

# resetting phase density    
with open("PR_delta=0_ell=4_kappa=1.csv", "r") as file:
   lines = csv.reader(file)
   for line in lines:
       x = float(line[0])
       P = float(line[1])
       reset_density_a.append(P)
       
# subfigure b

x_list_b = []
density_b = []
diff_density_b = []
reset_density_b = []

with open("P_delta=0.2_ell=4_kappa=1.csv", "r") as file:
   lines = csv.reader(file)
   for line in lines:
       x = float(line[0])
       P = float(line[1])
       x_list_b.append(x)
       density_b.append(P)
       
with open("PD_delta=0.2_ell=4_kappa=1.csv", "r") as file:
   lines = csv.reader(file)
   for line in lines:
       x = float(line[0])
       P = float(line[1])
       diff_density_b.append(P)
       
with open("PR_delta=0.2_ell=4_kappa=1.csv", "r") as file:
   lines = csv.reader(file)
   for line in lines:
       x = float(line[0])
       P = float(line[1])
       reset_density_b.append(P)
       
# subfigure c

x_list_c = []
density_c = []
diff_density_c = []
reset_density_c = []

with open("P_delta=0.5_ell=4_kappa=1.csv", "r") as file:
   lines = csv.reader(file)
   for line in lines:
       x = float(line[0])
       P = float(line[1])
       x_list_c.append(x)
       density_c.append(P)
       
with open("PD_delta=0.5_ell=4_kappa=1.csv", "r") as file:
   lines = csv.reader(file)
   for line in lines:
       x = float(line[0])
       P = float(line[1])
       diff_density_c.append(P)
       
with open("PR_delta=0.5_ell=4_kappa=1.csv", "r") as file:
   lines = csv.reader(file)
   for line in lines:
       x = float(line[0])
       P = float(line[1])
       reset_density_c.append(P)
       
# subfigure d

x_list_d = []
density_d = []
diff_density_d = []
reset_density_d = []

with open("P_delta=0_ell=1_kappa=4.csv", "r") as file:
   lines = csv.reader(file)
   for line in lines:
       x = float(line[0])
       P = float(line[1])
       x_list_d.append(x)
       density_d.append(P)
       
with open("PD_delta=0_ell=1_kappa=4.csv", "r") as file:
   lines = csv.reader(file)
   for line in lines:
       x = float(line[0])
       P = float(line[1])
       diff_density_d.append(P)
       
with open("PR_delta=0_ell=1_kappa=4.csv", "r") as file:
   lines = csv.reader(file)
   for line in lines:
       x = float(line[0])
       P = float(line[1])
       reset_density_d.append(P)

###############################################################################

############################## density figure #################################

###############################################################################

# define figure
fig1 = plt.figure(1, figsize=(16,16))
ax1 = fig1.add_subplot(221)

# theory plots (lines)
ax1.plot(x_list_a, density_a, color = 'black', linewidth = 5, label=r'$P(x)$')
ax1.plot(x_list_a, diff_density_a, color = 'blue', linewidth = 5, label=r'$P_{D}(x)$')
ax1.plot(x_list_a, reset_density_a, color = 'red', linewidth = 5, label=r'$P_{R}(x)$')

# axes labels
plt.xlabel(r'$x/L$', size = 25)
plt.ylabel(r'$P(x)$, $P_{D}(x)$, $P_{R}(x)$', size = 25, labelpad = 10)

# axes limits
plt.xlim(0, 1)
plt.ylim(0, )

# axes ticks
ax1.tick_params(which = 'both', direction='out', top=False, right=False)
ax1.tick_params(which = 'major', length = 6)
ax1.tick_params(which = 'minor', length = 3)

# axes tick locations
ax1.yaxis.set_minor_locator(MultipleLocator(0.05))
ax1.yaxis.set_major_locator(MultipleLocator(0.2))
ax1.xaxis.set_minor_locator(MultipleLocator(0.05))

# axes tick labels
plt.xticks(np.arange(0, 1.01, step=0.2), ('$0$', '$0.2$', '$0.4$', '$0.6$', '$0.8$', '$1$'), 
          size = 20)
plt.yticks(size=20)
ax1.tick_params(axis='x', which='major', pad=5)
ax1.tick_params(axis='y', which='major', pad=5)

# subfigure label
ax1.text(-0.21, 1.04, "(a)", size = 40, ha="left", va="top", transform=ax1.transAxes)

#legend
ax1.legend(loc = 'best', columnspacing = 0.7, ncol = 3, prop={'size': 20}, handletextpad=0.2,
           frameon=False)

# tight layout
plt.tight_layout()

plt.show()

###############################################################################

# define figure
ax2 = fig1.add_subplot(222)

# theory plots (lines)
ax2.plot(x_list_b, density_b, color = 'black', linewidth = 5, label=r'$P(x)$')
ax2.plot(x_list_b, diff_density_b, color = 'blue', linewidth = 5, label=r'$P_{D}(x)$')
ax2.plot(x_list_b, reset_density_b, color = 'red', linewidth = 5, label=r'$P_{R}(x)$')

# plot vertical line where ratchet peak appears
plt.vlines(0.2, 0, 1.4, linestyle = 'dashed', alpha = 0.1, linewidth = 3.5,
            color = 'black')

# axes labels
plt.xlabel(r'$x/L$', size = 25)
plt.ylabel(r'$P(x)$, $P_{D}(x)$, $P_{R}(x)$', size = 25, labelpad = 10)

# axes limits
plt.xlim(0, 1)
plt.ylim(0, 1.2)

# axes ticks
ax2.tick_params(which = 'both', direction='out', top=False, right=False)
ax2.tick_params(which = 'major', length = 6)
ax2.tick_params(which = 'minor', length = 3)

# axes tick locations
ax2.yaxis.set_minor_locator(MultipleLocator(0.05))
ax2.yaxis.set_major_locator(MultipleLocator(0.2))
ax2.xaxis.set_minor_locator(MultipleLocator(0.05))

# axes tick labels
plt.xticks(np.arange(0, 1.01, step=0.2), ('$0$', '$0.2$', '$0.4$', '$0.6$', '$0.8$', '$1$'), 
          size = 20)
plt.yticks(size=20)
ax2.tick_params(axis='x', which='major', pad=5)
ax2.tick_params(axis='y', which='major', pad=5)

# subfigure label
ax2.text(-0.21, 1.04, "(b)", size = 40, ha="left", va="top", transform=ax2.transAxes)

#legend
ax2.legend(loc = 'best', columnspacing = 0.7, ncol = 3, prop={'size': 20}, handletextpad=0.2,
           frameon=False)

# tight layout
plt.tight_layout()

plt.show()

###############################################################################

# define figure
ax3 = fig1.add_subplot(223)

# theory plots (lines)
ax3.plot(x_list_c, density_c, color = 'black', linewidth = 5, label=r'$P(x)$')
ax3.plot(x_list_c, diff_density_c, color = 'blue', linewidth = 5, label=r'$P_{D}(x)$')
ax3.plot(x_list_c, reset_density_c, color = 'red', linewidth = 5, label=r'$P_{R}(x)$')

# plot vertical line where ratchet peak appears
plt.vlines(0.5, 0, 1.4, linestyle = 'dashed', alpha = 0.1, linewidth = 3.5,
            color = 'black')

# axes labels
plt.xlabel(r'$x/L$', size = 25)
plt.ylabel(r'$P(x)$, $P_{D}(x)$, $P_{R}(x)$', size = 25, labelpad = 10)

# axes limits
plt.xlim(0, 1)
plt.ylim(0, 1.25)

# axes ticks
ax3.tick_params(which = 'both', direction='out', top=False, right=False)
ax3.tick_params(which = 'major', length = 6)
ax3.tick_params(which = 'minor', length = 3)

# axes tick locations
ax3.yaxis.set_minor_locator(MultipleLocator(0.05))
ax3.yaxis.set_major_locator(MultipleLocator(0.2))
ax3.xaxis.set_minor_locator(MultipleLocator(0.05))

# axes tick labels
plt.xticks(np.arange(0, 1.01, step=0.2), ('$0$', '$0.2$', '$0.4$', '$0.6$', '$0.8$', '$1$'), 
          size = 20)
plt.yticks(size=20)
ax3.tick_params(axis='x', which='major', pad=5)
ax3.tick_params(axis='y', which='major', pad=5)

# subfigure label
ax3.text(-0.21, 1.04, "(c)", size = 40, ha="left", va="top", transform=ax3.transAxes)

#legend
ax3.legend(loc = 'best', columnspacing = 0.7, ncol = 3, prop={'size': 20}, handletextpad=0.2,
           frameon=False)

# tight layout
plt.tight_layout()

plt.show()

###############################################################################

# define figure
ax4 = fig1.add_subplot(224)

# theory plots (lines)
ax4.plot(x_list_d, density_d, color = 'black', linewidth = 5, label=r'$P(x)$')
ax4.plot(x_list_d, diff_density_d, color = 'blue', linewidth = 5, label=r'$P_{D}(x)$')
ax4.plot(x_list_d, reset_density_d, color = 'red', linewidth = 5, label=r'$P_{R}(x)$')

# axes labels
plt.xlabel(r'$x/L$', size = 25)
plt.ylabel(r'$P(x)$, $P_{D}(x)$, $P_{R}(x)$', size = 25, labelpad = 10)

# axes limits
plt.xlim(0, 1)
plt.ylim(0, )

# axes ticks
ax4.tick_params(which = 'both', direction='out', top=False, right=False)
ax4.tick_params(which = 'major', length = 6)
ax4.tick_params(which = 'minor', length = 3)

# axes tick locations
ax4.yaxis.set_minor_locator(MultipleLocator(0.05))
ax4.yaxis.set_major_locator(MultipleLocator(0.2))
ax4.xaxis.set_minor_locator(MultipleLocator(0.05))

# axes tick labels
plt.xticks(np.arange(0, 1.01, step=0.2), ('$0$', '$0.2$', '$0.4$', '$0.6$', '$0.8$', '$1$'), 
          size = 20)
plt.yticks(size=20)
ax4.tick_params(axis='x', which='major', pad=5)
ax4.tick_params(axis='y', which='major', pad=5)

# subfigure label
ax4.text(-0.21, 1.04, "(d)", size = 40, ha="left", va="top", transform=ax4.transAxes)

#legend
ax4.legend(loc = 'center', columnspacing = 0.7, ncol = 3, prop={'size': 20}, handletextpad=0.2,
           frameon=False)

# tight layout
plt.tight_layout()

plt.show()

###############################################################################

end_time = time.time()

print('The full program took', end_time - start_time, 'seconds to run')