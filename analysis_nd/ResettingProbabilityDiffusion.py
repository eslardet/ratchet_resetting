import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
from density import *


small = 28
big = 30

plt.rc('font', size=40)          # controls default text sizes
plt.rc('axes', labelsize=big)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=small)    # fontsize of the tick labels
plt.rc('ytick', labelsize=small)    # fontsize of the tick labels
plt.rc('legend', fontsize=big)    # legend fontsize
markersize = 100

# matplotlib.rcParams["font.family"] = "serif"
plt.rcParams['text.usetex'] = True
# matplotlib.rcParams['mathtext.fontset'] = 'cm'
matplotlib.rcParams['font.family'] = 'STIXGeneral'



labels = ['(a)', '(b)', '(c)', '(d)']
alpha = 4.0
beta = 1.0
gamma = 0.5
bins = 20
t = format(20.0, '.6f')
t = format(0.0, '.6f')


phase = 'd'

# t_arr = t_arr = [format(i, '.6f') for i in np.arange(0, 1.01, 0.1)] + [format(20, '.6f')]

prob = get_prob_phase(alpha, beta, gamma, t, phase)

print(prob)


