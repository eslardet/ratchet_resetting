import os
import csv
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)


small = 28
big = 30

plt.rc('font', size=40)          # controls default text sizes
plt.rc('axes', labelsize=big)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=small)    # fontsize of the tick labels
plt.rc('ytick', labelsize=small)    # fontsize of the tick labels
plt.rc('legend', fontsize=big)    # legend fontsize

matplotlib.rcParams["font.family"] = 'STIXGeneral'
plt.rcParams['text.usetex'] = True


# colors = plt.cm.viridis(np.linspace(0, 1, 3))
# colors = plt.cm.Reds(np.linspace(0.2, 1, len(d_list)))

pi = np.pi
exp = np.exp
cos = np.cos
sin = np.sin
sqrt = np.sqrt
sinh = np.sinh
tanh = np.tanh
cosh = np.cosh


def P_D_Full(l,k,d):
    
    G_int = 2
    A = -k*(1 - exp(l))*(k + l*(d - 1))*(-d*l + k)*(d*l + k)*(exp(k/(1 - d)) - exp(k*(d - 2)/(d - 1)))*(-d*l + k + l)*exp(d*l + k/(d - 1))/(k*l*(-(-k + l*(d - 1))*(k + l*(2*d - 1))*(d*l - k)*exp(2*d*l) - (-k + l*(2*d - 1)*(k - 1))*(-d*l + k)*(-d*l + k + l)*exp(2*d*l + k) + (k + l*(d - 1))*(k + l*(2*d - 1)*(k - 1))*(d*l + k)*exp(k + l) - (k + l*(d - 1))*(d*l + k)*(-2*d*l + k + l)*exp(l)) + 2*(-k*(-d*l**4*(d - 1)*(k*(1 - 2*d)**2 - 1) - 2*k**4 + k**2*l**2*(2*d*(d - 1) + k*(1 - 2*d)**2 + 2) + (-d*l**4*(d - 1) + 2*k**4 - 2*k**2*l**2*(d*(d - 1) + 1))*cosh(k) + (d*l**4*(1 - 2*d)**2*(d - 1) + 2*k**4 - 2*k**2*l**2*(3*d*(d - 1) + 1))*sinh(k))*sinh(l/2) + l*(2*d**2*l**4*(d - 1)**2 + k**4 - k**2*l**2*(d*(d - 1) + k*(1 - 2*d)**2 + 1) + k**2*(-k**2 + l**2*(3*d*(d - 1) + 1))*sinh(k) - (2*d**2*l**4*(d - 1)**2 + k**4 + k**2*l**2*(-d**2 + d - 1))*cosh(k))*cosh(l/2))*exp(l*(d + 1/2)))
    
    return A*G_int

def energy_rate_analytic(l, k, d):
    # return (l*k*csch(l/2)*np.sinh(l*d/2)*np.sinh(l*(1-d)/2))/(d*(1-d))
    return (k*np.sinh(d*l/2)*np.sinh(l*(1-d)/2)) / (d*l*(1-d)*np.sinh(l/2)) * P_D_Full(l,k,d)

def get_sim_file(l, k, d, t):
    return './simulation_data/ell{}_kappa{}_delta{}/energy_{}'.format(l, k, d, t)





l = 4.0
k = 1.0
d = 0.2
t = format(100.0, '.6f')

file_name = get_sim_file(l, k, d, t)

with open(file_name) as f:
    reader = csv.reader(f, delimiter="\n")
    r = list(reader)
energy = np.array([float(i[0]) for i in r])

e_mean = np.round(np.mean(energy/float(t)),2)
e_st = np.round(np.std(energy/float(t)),2)


print("l = " + str(l) + ", k = " + str(k) + ", d = " + str(d) + ", t = " + str(t))
print("Simulation: " + str(e_mean) + " +- " + str(e_st))

print("Analytic: " + str(np.round(energy_rate_analytic(l, k, d),2)))