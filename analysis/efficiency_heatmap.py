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

matplotlib.rcParams["font.family"] = "serif"
plt.rcParams['text.usetex'] = True
# plt.rcParams['axes.labelpad']=10

l_list = np.linspace(0.1,10, 300)
k_list = np.linspace(0.01, 30, 300)

efficiency = np.zeros((len(l_list), len(k_list)))


fig, ax = plt.subplots(figsize=(10, 6))

def coth(x):
    return np.cosh(x)/np.sinh(x)

def eff(k, l):
    num = l**2 * k * (l * coth(l/2) - k*coth(k/2))**2
    denom = (l-k)*(l+k)*(2*(l-k)*(l+k) + l**3*coth(l/2)-l**2*k*coth(k/2))

    return num/denom

for i, l in enumerate(l_list):
    for j, k in enumerate(k_list):
        eta = eff(k, l)
        efficiency[i,j] = eta
    
for i, l in enumerate(l_list):
    max_val = np.max(efficiency[i,:])
    max_index = np.argmax(efficiency[i,:])
    # ax.plot(k_list[max_index], l, 'o', color='red', markersize=2)

pmesh = ax.pcolormesh(k_list, l_list, efficiency, shading='auto')
ax.set_xlabel(r'$\kappa$', labelpad=0)
ax.set_ylabel(r'$\ell$', labelpad=0)

ax.xaxis.set_minor_locator(MultipleLocator(2))
ax.yaxis.set_minor_locator(MultipleLocator(1))

cb = fig.colorbar(pmesh, ax=ax)
# cb.set_label(r'$\eta$')
ax.text(37,5,r'$\eta$')

folder = os.path.abspath('./plots/efficiency/')
if not os.path.exists(folder):
    os.makedirs(folder)
filename = os.path.join(folder, 'efficiency_meshplot300.pdf')
plt.savefig(filename, bbox_inches='tight')

# plt.show()
