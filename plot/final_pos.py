import os
import numpy as np
import matplotlib.pyplot as plt


def get_pos_file(h, a, D, r, t, phase):
    """
    Get path as string to position file for given parameters and time

    Input
    h : height of potential
    a : location of potential maximum
    D : diffusion coefficient
    r : reset rate
    t : time
    phase : 'd' or 'r'

    Returns
    sim_dir : path to position file
    """

    sim_dir = os.path.abspath('./simulation_data/h' + str(h) + '_a' + str(a) + '/D' + str(D) + '/r' + str(r)) 

    pos_file = os.path.join(sim_dir, phase + '_pos_' + str(t))
    return pos_file

def read_pos(file_name):
    """
    Read position file and return array of positions

    Input
    file_name : path to position file

    Returns
    x_all : array of positions

    """
    x_all = []
    file = open(file_name, 'r')
    for line in file.readlines():
        x = float(line)
        x_all.append(x)
    file.close()

    return np.array(x_all)

def plot_pos_phase(h, a, D, r, t, phase):
    """
    Plot position of particle over time

    Input
    x_all : array of positions
    h : height of potential
    a : location of potential maximum
    D : diffusion coefficient
    r : reset rate
    t : time
    phase : 'd' or 'r'

    Returns
    None
    """

    fig, ax = plt.subplots()
    pos_file = get_pos_file(h, a, D, r, t, phase)
    x_all = read_pos(pos_file)
    # ax.hist(x_all, bins=100, density=True)

    hist, bin_edges = np.histogram(x_all, bins=20, density=True)
    bin_centers = 0.5*(bin_edges[1:]+bin_edges[:-1])
    ax.plot(bin_centers, hist)

    ax.set_xlabel(r'$x$')
    ax.set_ylabel(r'$P_' + phase + r'(x)$')
    
    plt.show()

def plot_pos_total(h, a, D, r, t):
    """
    Plot position of particle over time

    Input
    x_all : array of positions
    h : height of potential
    a : location of potential maximum
    D : diffusion coefficient
    r : reset rate
    t : time

    Returns
    None
    """

    fig, ax = plt.subplots()
    pos_file_d = get_pos_file(h, a, D, r, t, 'd')
    x_all_d = read_pos(pos_file_d)
    pos_file_r = get_pos_file(h, a, D, r, t, 'r')
    x_all_r = read_pos(pos_file_r)
    x_all = np.concatenate((x_all_d, x_all_r))
    # ax.hist(x_all, bins=100, density=True)

    hist, bin_edges = np.histogram(x_all, bins=20, density=True)
    bin_centers = 0.5*(bin_edges[1:]+bin_edges[:-1])
    ax.plot(bin_centers, hist)

    ax.set_xlabel(r'$x$')
    ax.set_ylabel(r'$P(x)$')
    
    plt.show()

h = 1
a = 0.2
D = 0.1
r = 0.2
t = 20
phase = 'r'

plot_pos_phase(h, a, D, r, t, phase)
# plot_pos_total(h, a, D, r, t)

