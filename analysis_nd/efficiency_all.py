

import time
import datetime
start_time = time.time()
print('Program started at', datetime.datetime.now().time())

import os
from numba import jit
import numpy as np
import csv
import matplotlib.pyplot as plt
import matplotlib
import pickle
from matplotlib.lines import Line2D
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
from pylab import meshgrid,cm,imshow,contour,clabel,colorbar,axis,title,show
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from sympy.parsing.mathematica import parse_mathematica
from sympy import var
matplotlib.rcParams['mathtext.fontset'] = 'cm'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams['text.usetex'] = True
to_rgba = matplotlib.colors.to_rgba
pi = np.pi
exp = np.exp
cos = np.cos
sin = np.sin
sqrt = np.sqrt
sinh = np.sinh
cosh = np.cosh
tanh = np.tanh

###############################################################################

################################# functions ###################################

###############################################################################

### expression for current from Mathematica
# l,k,d = var('l k d')
# expr = "-(((1 - E^l) k^3 (E^(k + 2 l d) (k - l d) (l + k - l d) -        E^(l + k) (k + l (-1 + d)) (k + l d) -        E^(l d) (-1 +           2 d) (-((1 + E^l) l k) - (-1 + E^l) k^2 + (-1 + E^             l) l^2 (-1 + d) d)))/((-1 + E^       l) (l k (-E^l (k + l (-1 + d)) (l + k - 2 l d) (k + l d) -           E^(2 l d) (-k + l (-1 + d)) (-k + l d) (k + l (-1 + 2 d)) -           E^(k + 2 l d) (k - l d) (l + k - l d) (-k +              l (-1 + k) (-1 + 2 d)) +           E^(l + k) (k + l (-1 + d)) (k + l d) (k +              l (-1 + k) (-1 + 2 d))) +        2 E^(l (1/2 +            d)) (l Cosh[l/            2] (k^4 + 2 l^4 (-1 + d)^2 d^2 -              l^2 k^2 (1 + k (1 - 2 d)^2 + (-1 + d) d) - (k^4 +                 2 l^4 (-1 + d)^2 d^2 + l^2 k^2 (-1 + d - d^2)) Cosh[               k] + k^2 (-k^2 + l^2 (1 + 3 (-1 + d) d)) Sinh[k]) -           k Sinh[l/            2] (-2 k^4 - l^4 (-1 + k (1 - 2 d)^2) (-1 + d) d +              l^2 k^2 (2 + k (1 - 2 d)^2 + 2 (-1 + d) d) + (2 k^4 -                 l^4 (-1 + d) d - 2 l^2 k^2 (1 + (-1 + d) d)) Cosh[               k] + (2 k^4 + l^4 (1 - 2 d)^2 (-1 + d) d -                 2 l^2 k^2 (1 + 3 (-1 + d) d)) Sinh[k])))))"
# python_expr = parse_mathematica(expr)
# print(python_expr)

@jit(nopython=True)
def Current_Full(l,k,d):
    
    J = -k**3*(1 - exp(l))*(-(2*d - 1)*(d*l**2*(d - 1)*(exp(l) - 1) - k**2*(exp(l) - 1) - k*l*(exp(l) + 1))*exp(d*l) - (k + l*(d - 1))*(d*l + k)*exp(k + l) + (-d*l + k)*(-d*l + k + l)*exp(2*d*l + k))/((k*l*(-(-k + l*(d - 1))*(k + l*(2*d - 1))*(d*l - k)*exp(2*d*l) - (-k + l*(2*d - 1)*(k - 1))*(-d*l + k)*(-d*l + k + l)*exp(2*d*l + k) + (k + l*(d - 1))*(k + l*(2*d - 1)*(k - 1))*(d*l + k)*exp(k + l) - (k + l*(d - 1))*(d*l + k)*(-2*d*l + k + l)*exp(l)) + 2*(-k*(-d*l**4*(d - 1)*(k*(1 - 2*d)**2 - 1) - 2*k**4 + k**2*l**2*(2*d*(d - 1) + k*(1 - 2*d)**2 + 2) + (-d*l**4*(d - 1) + 2*k**4 - 2*k**2*l**2*(d*(d - 1) + 1))*cosh(k) + (d*l**4*(1 - 2*d)**2*(d - 1) + 2*k**4 - 2*k**2*l**2*(3*d*(d - 1) + 1))*sinh(k))*sinh(l/2) + l*(2*d**2*l**4*(d - 1)**2 + k**4 - k**2*l**2*(d*(d - 1) + k*(1 - 2*d)**2 + 1) + k**2*(-k**2 + l**2*(3*d*(d - 1) + 1))*sinh(k) - (2*d**2*l**4*(d - 1)**2 + k**4 + k**2*l**2*(-d**2 + d - 1))*cosh(k))*cosh(l/2))*exp(l*(d + 1/2)))*(exp(l) - 1))
    
    return J

@jit(nopython=True)
def P_D_Full(l,k,d):
    
    G_int = 2
    A = -k*(1 - exp(l))*(k + l*(d - 1))*(-d*l + k)*(d*l + k)*(exp(k/(1 - d)) - exp(k*(d - 2)/(d - 1)))*(-d*l + k + l)*exp(d*l + k/(d - 1))/(k*l*(-(-k + l*(d - 1))*(k + l*(2*d - 1))*(d*l - k)*exp(2*d*l) - (-k + l*(2*d - 1)*(k - 1))*(-d*l + k)*(-d*l + k + l)*exp(2*d*l + k) + (k + l*(d - 1))*(k + l*(2*d - 1)*(k - 1))*(d*l + k)*exp(k + l) - (k + l*(d - 1))*(d*l + k)*(-2*d*l + k + l)*exp(l)) + 2*(-k*(-d*l**4*(d - 1)*(k*(1 - 2*d)**2 - 1) - 2*k**4 + k**2*l**2*(2*d*(d - 1) + k*(1 - 2*d)**2 + 2) + (-d*l**4*(d - 1) + 2*k**4 - 2*k**2*l**2*(d*(d - 1) + 1))*cosh(k) + (d*l**4*(1 - 2*d)**2*(d - 1) + 2*k**4 - 2*k**2*l**2*(3*d*(d - 1) + 1))*sinh(k))*sinh(l/2) + l*(2*d**2*l**4*(d - 1)**2 + k**4 - k**2*l**2*(d*(d - 1) + k*(1 - 2*d)**2 + 1) + k**2*(-k**2 + l**2*(3*d*(d - 1) + 1))*sinh(k) - (2*d**2*l**4*(d - 1)**2 + k**4 + k**2*l**2*(-d**2 + d - 1))*cosh(k))*cosh(l/2))*exp(l*(d + 1/2)))
    
    return A*G_int

@jit(nopython=True)
def EnergyRate(l,k,d):
    
    Ei = (k*np.sinh(d*l/2)*np.sinh(l*(1-d)/2)) / (d*l*(1-d)*np.sinh(l/2)) * P_D_Full(l,k,d)

    return Ei

@jit(nopython=True)
def Efficiency(l,k,d):
    
    Ei = EnergyRate(l,k,d)
    J = Current_Full(l,k,d)
    
    eta = J**2 * l**2 / Ei

    if np.isnan(eta) or np.isinf(eta):
        return Efficiency(l,k,1-d)

    return eta

@jit(nopython=True)
def CurrentScript_KappaEll(d):
    
    # independent variables
    k_list = []
    l_list = []

    # dependent variable
    current_list = []

    # parameter ranges
    k_start = 0.01155
    k_end = 10
    k_step = 0.1

    l_start = 0.01
    l_end = 10
    l_step = 0.1
    
    for l in np.arange(l_start,l_end+l_step,l_step): # y values
        
        k_inner = []
        l_inner = []
        current_inner = []
        
        for k in np.arange(k_start,k_end+k_step,k_step): # x values
            
            # J = Current_Full(l,k,d)
            eta = Efficiency(l,k,d)
            
            k_inner.append(k)
            l_inner.append(l)
            current_inner.append(eta)
        
        #minimum_MSD_list.append(v_inner[MSD_inner.index(min(MSD_inner))])
        k_list.append(k_inner)
        l_list.append(l_inner)
        current_list.append(current_inner)
        
    return k_list, l_list, current_list

@jit(nopython=True)
def CurrentScript_DeltaEll(k):
    
    # independent variables
    d_list = []
    l_list = []

    # dependent variable
    current_list = []

    # parameter ranges
    d_start = 0.001155
    d_end = 1
    d_step = 0.02

    l_start = 0.01
    l_end = 10
    l_step = 0.1
    
    for l in np.arange(l_start,l_end+l_step,l_step): # y values
        
        d_inner = []
        l_inner = []
        current_inner = []
        
        for d in np.arange(d_start,d_end+d_step,d_step): # x values
            
            # J = Current_Full(l,k,d)
            eta = Efficiency(l,k,d)
            
            d_inner.append(d)
            l_inner.append(l)
            current_inner.append(eta)
        
        #minimum_MSD_list.append(v_inner[MSD_inner.index(min(MSD_inner))])
        d_list.append(d_inner)
        l_list.append(l_inner)
        current_list.append(current_inner)
        
    return d_list, l_list, current_list

@jit(nopython=True)
def CurrentScript_DeltaKappa(l):
    
    # independent variables
    d_list = []
    k_list = []

    # dependent variable
    current_list = []

    # parameter ranges
    d_start = 0.001155
    d_end = 1
    d_step = 0.02

    k_start = 0.01
    k_end = 10
    k_step = 0.1
    
    for k in np.arange(k_start,k_end+k_step,k_step): # y values
        
        d_inner = []
        k_inner = []
        current_inner = []
        
        for d in np.arange(d_start,d_end+d_step,d_step): # x values
            
            # J = Current_Full(l,k,d)
            eta = Efficiency(l,k,d)
            
            d_inner.append(d)
            k_inner.append(k)
            current_inner.append(eta)
        
        #minimum_MSD_list.append(v_inner[MSD_inner.index(min(MSD_inner))])
        d_list.append(d_inner)
        k_list.append(k_inner)
        current_list.append(current_inner)
        
    return d_list, k_list, current_list

###############################################################################

############################## run scripts ####################################

###############################################################################

k_list_kl, l_list_kl, J_list_kl = CurrentScript_KappaEll(0.2)
d_list_dl, l_list_dl, J_list_dl = CurrentScript_DeltaEll(4.001)
d_list_dk, k_list_dk, J_list_dk = CurrentScript_DeltaKappa(4.001)

#(a) kappa vs delta
max_eta = np.where(J_list_dk == np.nanmax(J_list_dk))
print('The maximum efficiency in (a) is', np.nanmax(J_list_dk), 'at kappa =', k_list_dk[max_eta[0][0]][max_eta[1][0]], 'and delta =', d_list_dk[max_eta[0][0]][max_eta[1][0]])

print(Efficiency(4, 8.4, 0.0001))

#(b) ell vs delta
max_eta = np.where(J_list_dl == np.nanmax(J_list_dl))
print('The maximum efficiency in (b) is', np.nanmax(J_list_dl), 'at ell =', l_list_dl[max_eta[0][0]][max_eta[1][0]], 'and delta =', d_list_dl[max_eta[0][0]][max_eta[1][0]])
print(Efficiency(6.4, 4, 0.0001))

#(c) ell vs kappa
max_eta = np.where(J_list_kl == np.nanmax(J_list_kl))
print('The maximum efficiency in (c) is', np.nanmax(J_list_kl), 'at kappa =', k_list_kl[max_eta[0][0]][max_eta[1][0]], 'and ell =', l_list_kl[max_eta[0][0]][max_eta[1][0]])
print(Efficiency(5.1, 6.7, 0.2))

###############################################################################

################################ figures ######################################

###############################################################################

# 2D colour plot

# define figure
fig1 = plt.figure(1, figsize=(20,5))

ax1 = fig1.add_subplot(131)

# plot contour
ratio_plot = ax1.pcolor(d_list_dk, k_list_dk, J_list_dk, cmap = 'viridis')# , norm=colors.LogNorm()

# colorbar
cbar = colorbar(ratio_plot, fraction = 0.045)
# setting fraction = 0.045 magically keeps the colorbar the same size as the plot
cbar.ax.tick_params(labelsize=20)
cbar.ax.set_title(r'$\eta$', size = 25, pad = 15)
    
# axes labels
plt.xlabel(r'$\delta$', size = 25, labelpad = 2)
plt.ylabel(r'$\kappa$', size = 25, labelpad = 2)

# axes ticks
ax1.tick_params(which = 'both', direction='out', top=False, right=False)
ax1.tick_params(which = 'major', length = 4)
ax1.tick_params(which = 'minor', length = 2)

# axes tick locations for log scale
ax1.xaxis.set_minor_locator(MultipleLocator(0.05))
ax1.yaxis.set_minor_locator(MultipleLocator(0.5))
ax1.xaxis.set_major_locator(MultipleLocator(0.2))
ax1.yaxis.set_major_locator(MultipleLocator(2))

# axes tick labels
plt.xticks(size = 20)
plt.yticks(size = 20)
ax1.tick_params(axis='x', which='major', pad=8)
ax1.tick_params(axis='y', which='major', pad=5)

# subfigure label
ax1.text(-0.25, 1.05, "(a)", size = 30, ha="left", va="top", transform=ax1.transAxes)

# plt.tight_layout()
#ax2.set_aspect('equal')
# plt.show()

###############################################################################

ax2 = fig1.add_subplot(132)

# plot contour
ratio_plot = ax2.pcolor(d_list_dl, l_list_dl, J_list_dl, cmap = 'viridis')# , norm=colors.LogNorm()

# colorbar
cbar = colorbar(ratio_plot, fraction = 0.045)
# setting fraction = 0.045 magically keeps the colorbar the same size as the plot
cbar.ax.tick_params(labelsize=20)
cbar.ax.set_title(r'$\eta$', size = 25, pad = 15)
    
# axes labels
plt.xlabel(r'$\delta$', size = 25, labelpad = 2)
plt.ylabel(r'$\ell$', size = 25, labelpad = 2)

# axes ticks
ax2.tick_params(which = 'both', direction='out', top=False, right=False)
ax2.tick_params(which = 'major', length = 4)
ax2.tick_params(which = 'minor', length = 2)

# axes tick locations for log scale
ax2.xaxis.set_minor_locator(MultipleLocator(0.05))
ax2.yaxis.set_minor_locator(MultipleLocator(0.5))
ax2.xaxis.set_major_locator(MultipleLocator(0.2))
ax2.yaxis.set_major_locator(MultipleLocator(2))

# axes tick labels
plt.xticks(size = 20)
plt.yticks(size = 20)
ax2.tick_params(axis='x', which='major', pad=8)
ax2.tick_params(axis='y', which='major', pad=5)

# subfigure label
ax2.text(-0.25, 1.05, "(b)", size = 30, ha="left", va="top", transform=ax2.transAxes)

# plt.tight_layout()
# plt.show()

###############################################################################

ax3 = fig1.add_subplot(133)

# plot contour
ratio_plot = ax3.pcolor(k_list_kl, l_list_kl, J_list_kl, cmap = 'viridis', vmin = 0)# , norm=colors.LogNorm()

# colorbar
cbar = colorbar(ratio_plot, fraction = 0.045)
# setting fraction = 0.045 magically keeps the colorbar the same size as the plot
cbar.ax.tick_params(labelsize=20)
cbar.ax.set_title(r'$\eta$', size = 25, pad = 15)
    
# axes labels
plt.xlabel(r'$\kappa$', size = 25, labelpad = 2)
plt.ylabel(r'$\ell$', size = 25, labelpad = 2)

# axes ticks
ax3.tick_params(which = 'both', direction='out', top=False, right=False)
ax3.tick_params(which = 'major', length = 4)
ax3.tick_params(which = 'minor', length = 2)

# axes tick locations for log scale
ax3.xaxis.set_minor_locator(MultipleLocator(0.5))
ax3.yaxis.set_minor_locator(MultipleLocator(0.5))
ax3.xaxis.set_major_locator(MultipleLocator(2))
ax3.yaxis.set_major_locator(MultipleLocator(2))

# axes tick labels
plt.xticks(size = 20)
plt.yticks(size = 20)
ax3.tick_params(axis='x', which='major', pad=8)
ax3.tick_params(axis='y', which='major', pad=5)

# subfigure label
ax3.text(-0.25, 1.05, "(c)", size = 30, ha="left", va="top", transform=ax3.transAxes)

plt.tight_layout(pad=2.0)
# plt.show()

folder = './plots/efficiency/'
file_name = 'EfficiencyVsParameters.pdf'

if not os.path.exists(folder):
    os.makedirs(folder)
plt.savefig(os.path.join(folder, file_name), bbox_inches='tight')
plt.close()

###############################################################################

end_time = time.time()

print('The full program took', end_time - start_time, 'seconds to run')