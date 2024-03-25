import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
from matplotlib import colors

small = 22
big = 28

plt.rc('font', size=big)          # controls default text sizes
plt.rc('axes', labelsize=big)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=small)    # fontsize of the tick labels
plt.rc('ytick', labelsize=small)    # fontsize of the tick labels
plt.rc('legend', fontsize=small)    # legend fontsize

# matplotlib.rcParams["font.family"] = "serif"
# plt.rcParams['text.usetex'] = True
matplotlib.rcParams['mathtext.fontset'] = 'cm'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
# plt.rcParams['axes.labelpad']=10


def coth(x):
    return np.cosh(x)/np.sinh(x)

def eff(k, l):
    num = k * (l * coth(l/2) - k*coth(k/2))**2
    denom = (l-k)*(l+k)*(2*(l-k)*(l+k) + l**3*coth(l/2)-l**2*k*coth(k/2))

    return num/denom

fig, axs = plt.subplots(1,2, figsize=(20, 6))


## (a) ##
ax = axs[0]
ax.text(-0.08, 0.95, r'(a)', horizontalalignment='center', verticalalignment='center', transform=ax.transAxes)

l_list = [0.1,1,2,5]
colors = plt.cm.Blues(np.linspace(0.2, 1, len(l_list)))

for i, l in enumerate(l_list):
    k= np.linspace(0.01, 30, 100)
    eta = eff(k, l)
    ax.plot(k, eta, label=r'$\ell={}$'.format(l), color=colors[i])

# l = 1
# k= np.linspace(0.01, 50, 100)
# eta = eff(k, l)
# ax.plot(k, eta, label=r'$\ell={}$'.format(l))

ax.set_xlabel(r'$\kappa$', labelpad=0)
ax.set_ylabel(r'$\eta$', labelpad=10)
ax.legend(frameon=False)

ax.set_xlim(0,30)
ax.xaxis.set_minor_locator(MultipleLocator(1))

ax.set_ylim(bottom=0)
ax.yaxis.set_minor_locator(MultipleLocator(0.002))



## (b) ##
ax = axs[1]
ax.text(-0.1, 0.95, r'(b)', horizontalalignment='center', verticalalignment='center', transform=ax.transAxes)

l_list = np.linspace(0.001,10, 500)
k_list = np.linspace(0.01, 30, 500)

efficiency = np.zeros((len(l_list), len(k_list)))

for i, l in enumerate(l_list):
    for j, k in enumerate(k_list):
        eta = eff(k, l)
        efficiency[i,j] = eta
    
for i, l in enumerate(l_list):
    max_val = np.max(efficiency[i,:])
    max_index = np.argmax(efficiency[i,:])
    ax.plot(k_list[max_index], l, 'o', color='red', markersize=2)

pmesh = ax.pcolormesh(k_list, l_list, efficiency, shading='auto', vmin=0)
ax.set_xlabel(r'$\kappa$', labelpad=0)
ax.set_ylabel(r'$\ell$', labelpad=0)


ax.xaxis.set_minor_locator(MultipleLocator(1))
ax.yaxis.set_minor_locator(MultipleLocator(1))
ax.tick_params(which = 'both', direction='out', top=False, right=False)
ax.tick_params(which = 'major', length = 6)
ax.tick_params(which = 'minor', length = 3)

ax.set_ylim(0, 10)

cbaxes = fig.add_axes([0.92, 0.12, 0.02, 0.75]) 
cb = fig.colorbar(pmesh, cax=cbaxes)
# cb.set_label(r'$\eta$')
ax.text(37,5,r'$\eta$')

folder = os.path.abspath('./plots/efficiency/')
if not os.path.exists(folder):
    os.makedirs(folder)
filename = os.path.join(folder, 'efficiency_both.pdf')
plt.savefig(filename, bbox_inches='tight')

# plt.show()