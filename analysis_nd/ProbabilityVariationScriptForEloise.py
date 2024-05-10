import time
import datetime
start_time = time.time()
print('Program started at', datetime.datetime.now().time())

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
tanh = np.tanh
cosh = np.cosh

###############################################################################

################################# load data ###################################

###############################################################################

d_list = []
diff_l4_k1 = []

with open("ProbDiff_ell=4_kappa=1.csv", "r") as file:
   lines = csv.reader(file)
   for line in lines:
       x = float(line[0])
       P = float(line[1])
       d_list.append(x)
       diff_l4_k1.append((P+1)/2)

diff_l4_k2 = []

with open("ProbDiff_ell=4_kappa=2.csv", "r") as file:
   lines = csv.reader(file)
   for line in lines:
       x = float(line[0])
       P = float(line[1])
       diff_l4_k2.append((P+1)/2)
       
diff_l5_k1 = []

with open("ProbDiff_ell=5_kappa=1.csv", "r") as file:
   lines = csv.reader(file)
   for line in lines:
       x = float(line[0])
       P = float(line[1])
       diff_l5_k1.append((P+1)/2)
       
diff_l5_k2 = []

with open("ProbDiff_ell=5_kappa=2.csv", "r") as file:
   lines = csv.reader(file)
   for line in lines:
       x = float(line[0])
       P = float(line[1])
       diff_l5_k2.append((P+1)/2)


###############################################################################

############################# Mathematica parser ##############################

###############################################################################

### F1 integral
# l,k,d = var('l k d')
# expr = "-((E^(k/(1 - d) - l  (1 + 2  d)) l  d  (E^(l  (2 + d)) k  (k + l  (-1 + d))  (l + k - 2  l  d)  (k + l  d) + E^(l + k +          2  l  d)  (k + l  (-1 + d))  (k - l  d)  (l + k - l  d)  (k +           l  d) +        E^(k + 2  l  (1 + d))  (k + l  (-1 + d))  (k - l  d)  (l + k -           l  d)  (k + l  d) +        E^(l + k + 3  l  d)          k  (k - l  d)  (l + k - l  d)  (l - k + l  (-2 + k)  d) -        E^(k + l  (2 + d))          k  (k + l  (-1 + d))  (k + l  d)  (l + k + l  (-2 + k)  d) +        E^(l - k +          2  l  d) (1 + E^l)  l^2  (-1 + d)  d  (k^2 +           l^2  (-1 + d)  d) +       E^(l + 3  l  d          k  (k - l  d)  (l + k - l  d)  (k + l  (-1 + 2  d))-       E^(2  l  (1 + d))  (-k + l  (-1 + d))  (k + l  d)  (-k^2 +           2  l^2  (-1 + d)  d +          l  k  (1 + d  (2 - k + 2  (-1 + k)  d))) -       E^(l + 2  l  d)  (k + l  (-1 + d))  (-k + l  d)  (-k^2 +           2  l^2  (-1 + d)  d +           l  k  (-1 + d  (-2 + k + 2  d - 2  k  d))) -        2  E^(-k + l  (3/2 + 2  d))          l  k  (-1 + d)  d  (-2  k^2 +           l^2  (1 + 2  (-1 + d)  d))  Sinh[l/2]))/((-1 + E^l)  (E^(k/(       1 - d)) - E^((k  (-2 + d))/(-1 + d)))  k  (k +        l  (-1 + d))  (k - l  d)  (l + k - l  d)  (k + l  d))))"
# python_expr = parse_mathematica(expr)
# Fone_int = -d*l*(d*l**2*(d - 1)*(d*l**2*(d - 1) + k**2)*(exp(l) + 1)*exp(2*d*l - k + l) + k*(k + l*(d - 1))*(d*l + k)*(-2*d*l + k + l)*exp(l*(d + 2)) - k*(k + l*(d - 1))*(d*l + k)*(d*l*(k - 2) + k + l)*exp(k + l*(d + 2)) + k*(-d*l + k)*(-d*l + k + l)*(d*l*(k - 2) - k + l)*exp(3*d*l + k + l) + (k + l*(d - 1))*(-d*l + k)*(d*l + k)*(-d*l + k + l)*exp(k + 2*l*(d + 1)) + (k + l*(d - 1))*(-d*l + k)*(d*l + k)*(-d*l + k + l)*exp(2*d*l + k + l) + exp(-2*d*k*l*(d - 1)*(-2*k**2 + l**2*(2*d*(d - 1) + 1))*exp(-k + l*(2*d + 3/2))*sinh(l/2) + 3*d*k*l*(k + l*(2*d - 1))*(-d*l + k)*(-d*l + k + l) + l - (-k + l*(d - 1))*(d*l + k)*(2*d*l**2*(d - 1) - k**2 + k*l*(d*(2*d*(k - 1) - k + 2) + 1))*exp(2*l*(d + 1)) - (k + l*(d - 1))*(d*l - k)*(2*d*l**2*(d - 1) - k**2 + k*l*(d*(-2*d*k + 2*d + k - 2) - 1))*exp(2*d*l + l)))*exp(k/(1 - d) - l*(2*d + 1))/(k*(k + l*(d - 1))*(-d*l + k)*(d*l + k)*(exp(l) - 1)*(exp(k/(1 - d)) - exp(k*(d - 2)/(d - 1)))*(-d*l + k + l))


### F2 integral
# l,k,d = var('l k d')
# expr = "(l (-1 +      d) (E^(-((k d)/(-1 + d)))       l (l + 2 k) (k + l (-1 + d)) (-1 + d) d (-k + l d) +      E^(l - (k d)/(-1 + d))       l (l - 2 k) (-k + l (-1 + d)) (-1 + d) d (k + l d) +      E^(l + k/(1 - d) - l d)       k (k + l (-1 + d)) (l + k - 2 l d) (k + l d) +      E^(l + k + k/(       1 - d)) (k + l (-1 + d)) (k - l d) (l + k - l d) (k + l d) +      E^((k (-2 + d))/(-1 +        d)) (k + l (-1 + d)) (k - l d) (l + k - l d) (k + l d) +      E^(k + k/(1 - d) + l d)       k (k - l d) (l + k - l d) (l - k - l k + l (-2 + k) d) -      E^(l + k + k/(1 - d) - l d)       k (k + l (-1 + d)) (k + l d) (l + k - l k + l (-2 + k) d) +      E^(k/(1 - d) + l d)       k (k - l d) (l + k - l d) (k + l (-1 + 2 d)) -      E^(k/(1 -        d)) (k + l (-1 + d)) (k - l d) (k^2 - 2 l^2 (-1 + d) d +         l k (1 + k + (2 - 3 k) d + 2 (-1 + k) d^2)) -      E^(l + k/(       1 - d)) (-k + l (-1 + d)) (k + l d) (-k^2 + 2 l^2 (-1 + d) d +         l k (1 + k + (2 - 3 k) d + 2 (-1 + k) d^2))))/((-1 + E^l) (E^(     k/(1 - d)) - E^((k (-2 + d))/(-1 + d))) k (k + l (-1 + d)) (k -      l d) (l + k - l d) (k + l d))"
# python_expr = parse_mathematica(expr)
# Ftwo_int = l*(d - 1)*(d*l*(d - 1)*(-2*k + l)*(-k + l*(d - 1))*(d*l + k)*exp(-d*k/(d - 1) + l) + d*l*(d - 1)*(k + l*(d - 1))*(2*k + l)*(d*l - k)*exp(-d*k/(d - 1)) + k*(k + l*(d - 1))*(d*l + k)*(-2*d*l + k + l)*exp(-d*l + k/(1 - d) + l) - k*(k + l*(d - 1))*(d*l + k)*(d*l*(k - 2) - k*l + k + l)*exp(-d*l + k + k/(1 - d) + l) + k*(k + l*(2*d - 1))*(-d*l + k)*(-d*l + k + l)*exp(d*l + k/(1 - d)) + k*(-d*l + k)*(-d*l + k + l)*(d*l*(k - 2) - k*l - k + l)*exp(d*l + k + k/(1 - d)) - (-k + l*(d - 1))*(d*l + k)*(2*d*l**2*(d - 1) - k**2 + k*l*(2*d**2*(k - 1) + d*(2 - 3*k) + k + 1))*exp(k/(1 - d) + l) + (k + l*(d - 1))*(-d*l + k)*(d*l + k)*(-d*l + k + l)*exp(k*(d - 2)/(d - 1)) + (k + l*(d - 1))*(-d*l + k)*(d*l + k)*(-d*l + k + l)*exp(k + k/(1 - d) + l) - (k + l*(d - 1))*(-d*l + k)*(-2*d*l**2*(d - 1) + k**2 + k*l*(2*d**2*(k - 1) + d*(2 - 3*k) + k + 1))*exp(k/(1 - d)))/(k*(k + l*(d - 1))*(-d*l + k)*(d*l + k)*(exp(l) - 1)*(exp(k/(1 - d)) - exp(k*(d - 2)/(d - 1)))*(-d*l + k + l))


### A normalisation
# l,k,d = var('l k d')
# expr = "-((E^(k/(-1 + d) +       l d) (1 - E^l) (E^(k/(1 - d)) - E^((k (-2 + d))/(-1 + d))) Sqrt[     r/D0] k (k + l (-1 + d)) (k - l d) (l + k - l d) (k +        l d))/(l k (-E^l (k + l (-1 + d)) (l + k - 2 l d) (k + l d) -         E^(2 l d) (-k + l (-1 + d)) (-k + l d) (k + l (-1 + 2 d)) -         E^(k + 2 l d) (k - l d) (l + k - l d) (-k +            l (-1 + k) (-1 + 2 d)) +         E^(l + k) (k + l (-1 + d)) (k + l d) (k +            l (-1 + k) (-1 + 2 d))) +      2 E^(l (1/2 +          d)) (l Cosh[l/          2] (k^4 + 2 l^4 (-1 + d)^2 d^2 -            l^2 k^2 (1 + k (1 - 2 d)^2 + (-1 + d) d) - (k^4 +               2 l^4 (-1 + d)^2 d^2 + l^2 k^2 (-1 + d - d^2)) Cosh[k] +            k^2 (-k^2 + l^2 (1 + 3 (-1 + d) d)) Sinh[k]) -        k Sinh[l/          2] (-2 k^4 - l^4 (-1 + k (1 - 2 d)^2) (-1 + d) d +            l^2 k^2 (2 + k (1 - 2 d)^2 + 2 (-1 + d) d) + (2 k^4 -               l^4 (-1 + d) d - 2 l^2 k^2 (1 + (-1 + d) d)) Cosh[             k] + (2 k^4 + l^4 (1 - 2 d)^2 (-1 + d) d -               2 l^2 k^2 (1 + 3 (-1 + d) d)) Sinh[k]))))"
# python_expr = parse_mathematica(expr)
# A = -k*(1 - exp(l))*(k + l*(d - 1))*(-d*l + k)*(d*l + k)*(exp(k/(1 - d)) - exp(k*(d - 2)/(d - 1)))*(-d*l + k + l)*exp(d*l + k/(d - 1))/(k*l*(-(-k + l*(d - 1))*(k + l*(2*d - 1))*(d*l - k)*exp(2*d*l) - (-k + l*(2*d - 1)*(k - 1))*(-d*l + k)*(-d*l + k + l)*exp(2*d*l + k) + (k + l*(d - 1))*(k + l*(2*d - 1)*(k - 1))*(d*l + k)*exp(k + l) - (k + l*(d - 1))*(d*l + k)*(-2*d*l + k + l)*exp(l)) + 2*(-k*(-d*l**4*(d - 1)*(k*(1 - 2*d)**2 - 1) - 2*k**4 + k**2*l**2*(2*d*(d - 1) + k*(1 - 2*d)**2 + 2) + (-d*l**4*(d - 1) + 2*k**4 - 2*k**2*l**2*(d*(d - 1) + 1))*cosh(k) + (d*l**4*(1 - 2*d)**2*(d - 1) + 2*k**4 - 2*k**2*l**2*(3*d*(d - 1) + 1))*sinh(k))*sinh(l/2) + l*(2*d**2*l**4*(d - 1)**2 + k**4 - k**2*l**2*(d*(d - 1) + k*(1 - 2*d)**2 + 1) + k**2*(-k**2 + l**2*(3*d*(d - 1) + 1))*sinh(k) - (2*d**2*l**4*(d - 1)**2 + k**4 + k**2*l**2*(-d**2 + d - 1))*cosh(k))*cosh(l/2))*exp(l*(d + 1/2)))

###############################################################################

################################### Scripts ###################################

###############################################################################

@jit(nopython=True)
def P_D_Full(l,k,d):
    
    G_int = 2
    A = -k*(1 - exp(l))*(k + l*(d - 1))*(-d*l + k)*(d*l + k)*(exp(k/(1 - d)) - exp(k*(d - 2)/(d - 1)))*(-d*l + k + l)*exp(d*l + k/(d - 1))/(k*l*(-(-k + l*(d - 1))*(k + l*(2*d - 1))*(d*l - k)*exp(2*d*l) - (-k + l*(2*d - 1)*(k - 1))*(-d*l + k)*(-d*l + k + l)*exp(2*d*l + k) + (k + l*(d - 1))*(k + l*(2*d - 1)*(k - 1))*(d*l + k)*exp(k + l) - (k + l*(d - 1))*(d*l + k)*(-2*d*l + k + l)*exp(l)) + 2*(-k*(-d*l**4*(d - 1)*(k*(1 - 2*d)**2 - 1) - 2*k**4 + k**2*l**2*(2*d*(d - 1) + k*(1 - 2*d)**2 + 2) + (-d*l**4*(d - 1) + 2*k**4 - 2*k**2*l**2*(d*(d - 1) + 1))*cosh(k) + (d*l**4*(1 - 2*d)**2*(d - 1) + 2*k**4 - 2*k**2*l**2*(3*d*(d - 1) + 1))*sinh(k))*sinh(l/2) + l*(2*d**2*l**4*(d - 1)**2 + k**4 - k**2*l**2*(d*(d - 1) + k*(1 - 2*d)**2 + 1) + k**2*(-k**2 + l**2*(3*d*(d - 1) + 1))*sinh(k) - (2*d**2*l**4*(d - 1)**2 + k**4 + k**2*l**2*(-d**2 + d - 1))*cosh(k))*cosh(l/2))*exp(l*(d + 1/2)))
    
    return A*G_int

@jit(nopython=True)
def P_D_Full_Script_kl(d):
    
    # independent variables
    k_list = []
    l_list = []

    # dependent variable
    PD_list = []

    # parameter ranges
    k_start = 0.0011
    k_end = 10
    k_step = 0.1

    l_start = 0.0009
    l_end = 20
    l_step = 0.2
    
    for l in np.arange(l_start,l_end+l_step,l_step): # y values
        
        k_inner = []
        l_inner = []
        PD_inner = []
        
        for k in np.arange(k_start,k_end+k_step,k_step): # x values
            
            PD = P_D_Full(l,k,d)
            
            k_inner.append(k)
            l_inner.append(l)
            PD_inner.append(PD)
        
        k_list.append(k_inner)
        l_list.append(l_inner)
        PD_list.append(PD_inner)
        
    return k_list, l_list, PD_list

@jit(nopython=True)
def P_D_Full_Script_dl(k):
    
    # independent variables
    d_list = []
    l_list = []

    # dependent variable
    PD_list = []

    # parameter ranges
    d_start = 0.001155
    d_end = 1
    d_step = 0.0093

    l_start = 0.0009
    l_end = 10
    l_step = 0.1
    
    for l in np.arange(l_start,l_end+l_step,l_step): # y values
        
        d_inner = []
        l_inner = []
        PD_inner = []
        
        for d in np.arange(d_start,d_end+d_step,d_step): # x values
            
            PD = P_D_Full(l,k,d)
            if np.isnan(PD) == True:
                PD = P_D_Full(l,k,1-d)
            
            d_inner.append(d)
            l_inner.append(l)
            PD_inner.append(PD)
        
        d_list.append(d_inner)
        l_list.append(l_inner)
        PD_list.append(PD_inner)
        
    return d_list, l_list, PD_list

@jit(nopython=True)
def P_D_Full_Script_dk(l):
    
    # independent variables
    d_list = []
    k_list = []

    # dependent variable
    PD_list = []

    # parameter ranges
    d_start = 0.001155
    d_end = 1
    d_step = 0.0093

    k_start = 0.0009
    k_end = 10
    k_step = 0.1
    
    for k in np.arange(k_start,k_end+k_step,k_step): # y values
        
        d_inner = []
        k_inner = []
        PD_inner = []
        
        for d in np.arange(d_start,d_end+d_step,d_step): # x values
            
            PD = P_D_Full(l,k,d)
            if np.isnan(PD) == True or np.isinf(PD) == True:
                PD = P_D_Full(l,k,1-d)
            
            d_inner.append(d)
            k_inner.append(k)
            PD_inner.append(PD)
        
        d_list.append(d_inner)
        k_list.append(k_inner)
        PD_list.append(PD_inner)
        
    return d_list, k_list, PD_list

###############################################################################

############################## run scripts ####################################

###############################################################################

k_list_kl, l_list_kl, PD_list_kl = P_D_Full_Script_kl(0.2)
d_list_dl, l_list_dl, PD_list_dl = P_D_Full_Script_dl(4)
d_list_dk, k_list_dk, PD_list_dk = P_D_Full_Script_dk(8)
print(PD_list_dl[-1][-1])

###############################################################################

################################ figures ######################################

###############################################################################

# # 2D colour plot

# define figure
fig2 = plt.figure(2, figsize=(20,5))

ax1 = fig2.add_subplot(131)

# plot contour
dk_plot = ax1.pcolor(d_list_dk, k_list_dk, PD_list_dk, cmap = 'bwr_r', vmin = 0, vmax = 1)# , norm=colors.LogNorm()

# colorbar
cbar = colorbar(dk_plot, fraction = 0.045, ticks = [0,0.25,0.5,0.75,1])
#cbar.set_label(r'$p_{D} - p_{R}$', rotation = 0, size = 20, labelpad = 35)
# setting fraction = 0.045 magically keeps the colorbar the same size as the plot
cbar.ax.tick_params(labelsize=20)
cbar.ax.set_title(r'$p_{D}$', size = 25, pad = 25)
    
# axes labels
plt.xlabel(r'$\delta$', size = 25, labelpad = 2)
plt.ylabel(r'$\kappa$', size = 25, labelpad = 2)

# axes limits
plt.xlim(0, 1)

# axes ticks
ax1.tick_params(which = 'both', direction='out', top=False, right=False)
ax1.tick_params(which = 'major', length = 4)
ax1.tick_params(which = 'minor', length = 2)

# axes tick locations for log scale
#ax2.xaxis.set_minor_locator(MultipleLocator(0.05))
ax1.yaxis.set_minor_locator(MultipleLocator(0.5))
#ax2.xaxis.set_major_locator(MultipleLocator(0.2))
ax1.yaxis.set_major_locator(MultipleLocator(2))

# axes tick labels
plt.xticks(size = 20)
plt.yticks(size = 20)
ax1.tick_params(axis='x', which='major', pad=8)
ax1.tick_params(axis='y', which='major', pad=5)

# subfigure label
ax1.text(-0.25, 1.05, "(a)", size = 30, ha="left", va="top", transform=ax1.transAxes)

plt.tight_layout()
plt.show()

###############################################################################

ax2 = fig2.add_subplot(132)

# plot contour
dl_plot = ax2.pcolor(d_list_dl, l_list_dl, PD_list_dl, cmap = 'bwr_r', vmin = 0, vmax = 1)# , norm=colors.LogNorm()

# colorbar
cbar = colorbar(dl_plot, fraction = 0.045, ticks = [0,0.25,0.5,0.75,1])
#cbar.set_label(r'$p_{D} - p_{R}$', rotation = 0, size = 20, labelpad = 35)
# setting fraction = 0.045 magically keeps the colorbar the same size as the plot
cbar.ax.tick_params(labelsize=20)
cbar.ax.set_title(r'$p_{D}$', size = 25, pad = 25)
    
# axes labels
plt.xlabel(r'$\delta$', size = 25, labelpad = 2)
plt.ylabel(r'$\ell$', size = 25, labelpad = 2)

# axes limits
plt.xlim(0, 1)

# axes ticks
ax2.tick_params(which = 'both', direction='out', top=False, right=False)
ax2.tick_params(which = 'major', length = 4)
ax2.tick_params(which = 'minor', length = 2)

# axes tick locations for log scale
#ax2.xaxis.set_minor_locator(MultipleLocator(0.05))
ax2.yaxis.set_minor_locator(MultipleLocator(0.5))
#ax2.xaxis.set_major_locator(MultipleLocator(0.2))
ax2.yaxis.set_major_locator(MultipleLocator(2))

# axes tick labels
plt.xticks(size = 20)
plt.yticks(size = 20)
ax2.tick_params(axis='x', which='major', pad=8)
ax2.tick_params(axis='y', which='major', pad=5)

# subfigure label
ax2.text(-0.25, 1.05, "(b)", size = 30, ha="left", va="top", transform=ax2.transAxes)

plt.tight_layout()
plt.show()

###############################################################################

ax3 = fig2.add_subplot(133)

# plot contour
kl_plot = ax3.pcolor(k_list_kl, l_list_kl, PD_list_kl, cmap = 'bwr_r', vmin = 0, vmax = 1)# , norm=colors.LogNorm()

# colorbar
cbar = colorbar(kl_plot, fraction = 0.045, ticks = [0,0.25,0.5,0.75,1])
#cbar.set_label(r'$p_{D} - p_{R}$', rotation = 0, size = 20, labelpad = 35)
# setting fraction = 0.045 magically keeps the colorbar the same size as the plot
cbar.ax.tick_params(labelsize=20)
cbar.ax.set_title(r'$p_{D}$', size = 25, pad = 25)
    
# axes labels
plt.xlabel(r'$\kappa$', size = 25, labelpad = 2)
plt.ylabel(r'$\ell$', size = 25, labelpad = 2)

# axes limits
#plt.xlim(0, N)
#plt.ylim(0, M)

# axes ticks
ax3.tick_params(which = 'both', direction='out', top=False, right=False)
ax3.tick_params(which = 'major', length = 4)
ax3.tick_params(which = 'minor', length = 2)

# axes tick locations for log scale
ax3.xaxis.set_minor_locator(MultipleLocator(0.5))
ax3.yaxis.set_minor_locator(MultipleLocator(1))
ax3.xaxis.set_major_locator(MultipleLocator(2))
ax3.yaxis.set_major_locator(MultipleLocator(4))

# axes tick labels
plt.xticks(size = 20)
plt.yticks(size = 20)
ax3.tick_params(axis='x', which='major', pad=8)
ax3.tick_params(axis='y', which='major', pad=5)

# subfigure label
ax3.text(-0.25, 1.05, "(c)", size = 30, ha="left", va="top", transform=ax3.transAxes)

plt.tight_layout()
plt.show()

###############################################################################

end_time = time.time()

print('The full program took', end_time - start_time, 'seconds to run')