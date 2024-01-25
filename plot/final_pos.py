import os
import numpy as np
import matplotlib.pyplot as plt


def get_pos_file(h, a, D, r, t):
    """
    Get path as string to position file for given parameters and time

    Input
    h : height of potential
    a : location of potential maximum
    D : diffusion coefficient
    r : reset rate
    t : time

    Returns
    sim_dir : path to position file
    """

    sim_dir = os.path.abspath('./simulation_data/h' + str(h) + '_a' + str(a) + '/D' + str(D) + '/r' + str(r)) 

    pos_file = os.path.join(sim_dir, 'pos_' + str(t))
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

def plot_pos(h, a, D, r, t):
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
    pos_file = get_pos_file(h, a, D, r, t)
    x_all = read_pos(pos_file)
    ax.hist(x_all, bins=100)
    
    plt.show()

h = 1
a = 0.2
D = 0.1
r = 0.2
t = 10

plot_pos(h, a, D, r, t)