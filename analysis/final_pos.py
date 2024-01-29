import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib


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

def initalize_plotting(cols='discrete', num=10):
    small = 14
    big = 20

    plt.rc('font', size=small)          # controls default text sizes
    plt.rc('axes', labelsize=big)     # fontsize of the x and y labels
    plt.rc('xtick', labelsize=big)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=big)    # fontsize of the tick labels
    plt.rc('legend', fontsize=small)    # legend fontsize

    matplotlib.rcParams["font.family"] = "serif"
    plt.rcParams['text.usetex'] = True

    fig, ax = plt.subplots(figsize=(7,5))

    if cols == 'discrete':
        ax.set_prop_cycle(color=plt.cm.tab10.colors)
    elif cols == 'continuous':
        ax.set_prop_cycle(color=plt.cm.viridis(np.linspace(0, 1, num)))
    return fig, ax

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

def plot_pos_phase(h, a_arr, D, r_arr, t_arr, phase, bins=20, save_plot=False, show_plot=True):
    """
    Plot position of particle over time

    Input
    x_all : array of positions
    h : height of potential
    a_arr : locations of potential maximum
    D : diffusion coefficient
    r_arr : reset rates
    t_arr : times to plot
    phase : 'd' or 'r'

    Returns
    None
    """
    fig, ax = initalize_plotting(cols='continuous', num=len(a_arr)*len(r_arr)*len(t_arr))
    for a in a_arr:
        for r in r_arr:
            for t in t_arr:
                pos_file = get_pos_file(h, a, D, r, t, phase)
                x_all = read_pos(pos_file)
                # ax.hist(x_all, bins=100, density=True)

                hist, bin_edges = np.histogram(x_all, bins=bins, density=True)
                bin_centers = 0.5*(bin_edges[1:]+bin_edges[:-1])
                ax.plot(bin_centers, hist, label=r'$a={a},\  r={r},\ t={t}$'.format(a=a, r=r, t=t))

    ax.set_xlim(0,1)
    ax.set_xlabel(r'$x$')
    ax.set_ylabel(r'$P_' + phase + r'(x)$')
    ax.legend()
    
    if save_plot:
        folder = os.path.abspath('./plots/density/')
        if not os.path.exists(folder):
            os.makedirs(folder)
        file_name = os.path.join(folder, phase + '_h' + str(h) + '_D' + str(D) + '_r' + str(r) + '.pdf')
        plt.savefig(file_name, bbox_inches='tight')
        plt.close()
    if show_plot:
        plt.show()

def plot_pos_both(h, a, D, r, t, bins=20, save_plot=False, show_plot=True):
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

    fig, ax = initalize_plotting(cols='discrete')
    for phase in ['d', 'r']:
        pos_file = get_pos_file(h, a, D, r, t, phase)
        x_all = read_pos(pos_file)
        # ax.hist(x_all, bins=100, density=True)
        hist, bin_edges = np.histogram(x_all, bins=bins, density=True)
        bin_centers = 0.5*(bin_edges[1:]+bin_edges[:-1])
        ax.plot(bin_centers, hist, label="Phase " + phase)

    ax.set_title(r'$a={a}, \ r={r}, \ t={t}$'.format(a=a, r=r, t=t))
    ax.set_xlim(0,1)
    ax.set_xlabel(r'$x$')
    ax.set_ylabel(r'$P(x)$')
    ax.legend()
    
    if save_plot:
        folder = os.path.abspath('./plots/density/')
        if not os.path.exists(folder):
            os.makedirs(folder)
        file_name = os.path.join(folder, 'both_h' + str(h) + '_a' + str(a) + '_D' + str(D) + '_r' + str(r) + '.pdf')
        plt.savefig(file_name, bbox_inches='tight')
        plt.close()
    if show_plot:
        plt.show()

def plot_pos_total(h, a_arr, D, r_arr, t_arr, bins=20, save_plot=False, show_plot=True):
    """
    Plot position of particle over time

    Input
    x_all : array of positions
    h : height of potential
    a_arr : locations of potential maximum
    D : diffusion coefficient
    r_arr : reset rates
    t : times to plot

    Returns
    None
    """

    fig, ax = initalize_plotting(cols='continuous', num=len(a_arr)*len(r_arr)*len(t_arr))
    for a in a_arr:
        for r in r_arr:
            for t in t_arr:
                pos_file_d = get_pos_file(h, a, D, r, t, 'd')
                x_all_d = read_pos(pos_file_d)
                pos_file_r = get_pos_file(h, a, D, r, t, 'r')
                x_all_r = read_pos(pos_file_r)
                x_all = np.concatenate((x_all_d, x_all_r))
                # ax.hist(x_all, bins=100, density=True)

                hist, bin_edges = np.histogram(x_all, bins=bins, density=True)
                bin_centers = 0.5*(bin_edges[1:]+bin_edges[:-1])
                ax.plot(bin_centers, hist, label=r'$a={a}, \ r={r}, \ t={t}$'.format(a=a, r=r, t=t))

    ax.set_xlim(0,1)
    ax.set_xlabel(r'$x$')
    ax.set_ylabel(r'$P(x)$')
    ax.legend()
    
    if save_plot:
        folder = os.path.abspath('./plots/density/')
        if not os.path.exists(folder):
            os.makedirs(folder)
        file_name = os.path.join(folder, 'joint_h' + str(h) + '_D' + str(D) + '_r' + str(r) + '.pdf')
        plt.savefig(file_name, bbox_inches='tight')
        plt.close()
    if show_plot:
        plt.show()

h = 1
a = 0.2
a_arr = np.round(np.arange(0.1,0.6,0.1),1)
D = 0.1
r = 0.2
r_arr = [0.2]
t = 20
t_arr = np.arange(20, 21, 1)
phase = 'r'

plot_pos_phase(h, a_arr, D, r_arr, t_arr, phase, save_plot=True, show_plot=False)
# plot_pos_both(h, a, D, r, t, save_plot=True, show_plot=False)
# plot_pos_total(h, a, D, r, t_arr)

