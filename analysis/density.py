import os
import csv
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

def get_prob_phase(h, a, D, r, t, phase='r'):
    pos_file_d = get_pos_file(h, a, D, r, t, phase='d')
    x_all_d = read_pos(pos_file_d)
    pos_file_r = get_pos_file(h, a, D, r, t, phase='r')
    x_all_r = read_pos(pos_file_r)

    prob_d = len(x_all_d)/(len(x_all_d) + len(x_all_r))
    prob_r = len(x_all_r)/(len(x_all_d) + len(x_all_r))

    if phase == 'd':
        return prob_d
    elif phase == 'r':
        return prob_r

def plot_pos_phase(h, a_arr, D, r_arr, t_arr, phase, bins=20, save_plot=False, show_plot=True):
    """
    Plot position of particle over time/ array of a or r values for a certain phase
    """
    fig, ax = initalize_plotting(cols='continuous', num=len(a_arr)*len(r_arr)*len(t_arr))
    for a in a_arr:
        for r in r_arr:
            for t in t_arr:
                pos_file = get_pos_file(h, a, D, r, t, phase)
                x_all = read_pos(pos_file)

                prob = get_prob_phase(h, a, D, r, t, phase)
                # ax.hist(x_all, bins=100, density=True)
                hist, bin_edges = np.histogram(x_all, bins=bins, density=True)
                bin_centers = 0.5*(bin_edges[1:]+bin_edges[:-1])
                ax.plot(bin_centers, prob*hist, label=r'$a={a},\  r={r},\ t={t}$'.format(a=a, r=r, t=str(float(t))))

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
    Plot position of particle over time for each phase
    """

    fig, ax = initalize_plotting(cols='discrete')

    for phase in ['d', 'r']:
        pos_file = get_pos_file(h, a, D, r, t, phase)
        x_all = read_pos(pos_file)
        
        prob = get_prob_phase(h, a, D, r, t, phase)
        # ax.hist(x_all, bins=100, density=True)
        hist, bin_edges = np.histogram(x_all, bins=bins, density=True)
        bin_centers = 0.5*(bin_edges[1:]+bin_edges[:-1])
        ax.plot(bin_centers, prob*hist, label="Phase " + phase)


    ax.set_title(r'$a={a}, \ r={r}, \ t={t}$'.format(a=a, r=r, t=str(float(t))))
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
    Plot position of particle over time (combined phases)
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
                ax.plot(bin_centers, hist, label=r'$a={a}, \ r={r}, \ t={t}$'.format(a=a, r=r, t=str(float(t))))

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



def read_csv(filename):
    """
    Read csv output for plotting
    """
    folder = os.path.abspath('./analytical_data/')
    file = os.path.join(folder, filename)
    x_plot = []
    y_plot = []
    with open(file) as f:
        reader = csv.reader(f, delimiter=",")
        for line in reader:
            x_plot.append(float(line[0]))
            y_plot.append(float(line[1]))
    
    return np.array(x_plot), np.array(y_plot)

def plot_csv(filename, ax=None):
    """
    Plot csv output
    """
    x_plot, y_plot = read_csv(filename)
    if ax is None:
        fig, ax = initalize_plotting()
    ax.plot(x_plot, y_plot)
    # ax.set_xlim(0,1)
    # ax.set_xlabel(r'$x$')
    # ax.set_ylabel(r'$P(x)$')
    
    return ax

def plot_compare(h, a, D, r, t, phase, bins=20, ax=None):
    """
    Plot analytical and simulation data
    """
    filename = 'P_{}_h{}_a{}_D{}_r{}.csv'.format(phase, h, a, D, r)
    x_plot, y_plot = read_csv(filename)
    if ax is None:
        fig, ax = initalize_plotting()
    ax.plot(x_plot, y_plot, label='Analytical')

    pos_file = get_pos_file(h, a, D, r, t, phase)
    x_all = read_pos(pos_file)
    prob = get_prob_phase(h, a, D, r, t, phase)
    
    hist, bin_edges = np.histogram(x_all, bins=bins, density=True)
    bin_centers = 0.5*(bin_edges[1:]+bin_edges[:-1])
    ax.plot(bin_centers, prob*hist, 'o', label='Simulation')
    ax.set_xlim(0,1)
    ax.set_xlabel(r'$x$')
    ax.set_ylabel(r'$P(x)$')
    ax.legend()
    
    return ax