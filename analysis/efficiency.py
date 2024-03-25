import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)

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

l_list = [0.0]

# colors = plt.cm.viridis(np.linspace(0, 1, 3))
colors = plt.cm.Blues(np.linspace(0.2, 1, len(l_list)))



fig, ax = plt.subplots(figsize=(10, 6))

def coth(x):
    return np.cosh(x)/np.sinh(x)

def eff(k, l):
    num = k * (l * coth(l/2) - k*coth(k/2))**2
    denom = (l-k)*(l+k)*(2*(l-k)*(l+k) + l**3*coth(l/2)-l**2*k*coth(k/2))

    return num/denom

for i, l in enumerate(l_list):
    k= np.linspace(0.01, 50, 100)
    eta = eff(k, l)
    ax.plot(k, eta, label=r'$\ell={}$'.format(l), color=colors[i])

# l = 1
# k= np.linspace(0.01, 50, 100)
# eta = eff(k, l)
# ax.plot(k, eta, label=r'$\ell={}$'.format(l))

ax.set_xlabel(r'$\kappa$', labelpad=0)
ax.set_ylabel(r'$\eta$', labelpad=10)
ax.legend(frameon=False)

ax.set_xlim(0,50)
ax.xaxis.set_minor_locator(MultipleLocator(2))

ax.set_ylim(bottom=0)
ax.yaxis.set_minor_locator(MultipleLocator(0.002))


folder = os.path.abspath('./plots/efficiency/')
if not os.path.exists(folder):
    os.makedirs(folder)
filename = os.path.join(folder, 'efficiency_various_l.pdf')
# plt.savefig(filename, bbox_inches='tight')

plt.show()
