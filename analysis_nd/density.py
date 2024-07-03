import os
import csv
import numpy as np
import matplotlib.pyplot as plt
import matplotlib


def get_pos_file(ell, kappa, delta, t, phase, model='v1'):
    """
    Get path as string to position file for given parameters and time

    Input
    ell : ell value
    kappa : kappa value
    delta : delta value
    t : time
    phase : 'd' or 'r'

    Returns
    sim_dir : path to position file
    """

    if model == 'v1':
        sim_dir = os.path.abspath('./simulation_data/alpha{}_beta{}_gamma{}'.format(ell, kappa, delta)) 
    elif model == 'v2':
        sim_dir = os.path.abspath('./simulation_data_stop_start/alpha{}_beta{}_gamma{}'.format(ell, kappa, delta))

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

    matplotlib.rcParams["font.family"] = 'STIXGeneral'
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

def get_prob_phase(ell, kappa, delta, t, phase='r', model="v1"):
    pos_file_d = get_pos_file(ell, kappa, delta, t, phase='d', model=model)
    x_all_d = read_pos(pos_file_d)
    pos_file_r = get_pos_file(ell, kappa, delta, t, phase='r', model=model)
    x_all_r = read_pos(pos_file_r)

    prob_d = len(x_all_d)/(len(x_all_d) + len(x_all_r))
    prob_r = len(x_all_r)/(len(x_all_d) + len(x_all_r))

    if phase == 'd':
        return prob_d
    elif phase == 'r':
        return prob_r

def plot_pos_phase(ell, kappa_arr, delta_arr, t_arr, phase, bins=20, save_plot=False, show_plot=True):
    """
    Plot position of particle over time/ array of a or r values for a certain phase
    """
    fig, ax = initalize_plotting(cols='continuous', num=len(kappa_arr)*len(delta_arr)*len(t_arr))
    for kappa in kappa_arr:
        for delta in delta_arr:
            for t in t_arr:
                pos_file = get_pos_file(ell, kappa, delta, t, phase)
                x_all = read_pos(pos_file)

                prob = get_prob_phase(ell, kappa, delta, t, phase)
                # ax.hist(x_all, bins=100, density=True)
                hist, bin_edges = np.histogram(x_all, bins=bins, density=True)
                bin_centers = 0.5*(bin_edges[1:]+bin_edges[:-1])
                ax.plot(bin_centers, prob*hist, label=r'$\kappa={kappa},\  \delta={delta},\ t={t}$'.format(kappa=kappa, delta=delta, t=str(float(t))))

    ax.set_xlim(0,ell)
    ax.set_xlabel(r'$x$')
    ax.set_ylabel(r'$P_' + phase + r'(x)$')
    ax.legend()
    
    if save_plot:
        folder = os.path.abspath('./plots/density/')
        if not os.path.exists(folder):
            os.makedirs(folder)
        file_name = os.path.join(folder, phase + '_ell{}_kappa{}_delta{}.pdf'.format(ell, kappa, delta))
        plt.savefig(file_name, bbox_inches='tight')
        plt.close()
    if show_plot:
        plt.show()

def plot_pos_both(ell, kappa, delta, t, bins=20, save_plot=False, show_plot=True):
    """
    Plot position of particle over time for each phase
    """

    fig, ax = initalize_plotting(cols='discrete')

    for phase in ['d', 'r']:
        pos_file = get_pos_file(ell, kappa, delta, t, phase)
        x_all = read_pos(pos_file)
        
        prob = get_prob_phase(ell, kappa, delta, t, phase)
        # ax.hist(x_all, bins=100, density=True)
        hist, bin_edges = np.histogram(x_all, bins=bins, density=True)
        bin_centers = 0.5*(bin_edges[1:]+bin_edges[:-1])
        ax.plot(bin_centers, prob*hist, label="Phase " + phase)


    ax.set_title(r'$\ell={ell}, \ \kappa={kappa}, \delta={delta} \ t={t}$'.format(ell=ell, kappa=kappa, delta=delta, t=str(float(t))))
    ax.set_xlim(0,1)
    ax.set_xlabel(r'$x$')
    ax.set_ylabel(r'$P(x)$')
    ax.legend()
    
    if save_plot:
        folder = os.path.abspath('./plots/density/')
        if not os.path.exists(folder):
            os.makedirs(folder)
        file_name = os.path.join(folder, 'both_ell{}_kappa{}_delta{}_t{}.pdf'.format(ell, kappa, delta, t))
        plt.savefig(file_name, bbox_inches='tight')
        plt.close()
    if show_plot:
        plt.show()

def plot_pos_total(ell, kappa_arr, delta_arr, t_arr, bins=20, save_plot=False, show_plot=True):
    """
    Plot position of particle over time (combined phases)
    """

    fig, ax = initalize_plotting(cols='continuous', num=len(kappa_arr)*len(delta_arr)*len(t_arr))
    for kappa in kappa_arr:
        for delta in delta_arr:
            for t in t_arr:
                pos_file_d = get_pos_file(ell, kappa, delta, t, 'd')
                x_all_d = read_pos(pos_file_d)
                pos_file_r = get_pos_file(ell, kappa, delta, t, 'r')
                x_all_r = read_pos(pos_file_r)
                x_all = np.concatenate((x_all_d, x_all_r))
                # ax.hist(x_all, bins=100, density=True)

                hist, bin_edges = np.histogram(x_all, bins=bins, density=True)
                bin_centers = 0.5*(bin_edges[1:]+bin_edges[:-1])
                ax.plot(bin_centers, hist, label=r'$\kappa=$' + str(kappa) + r'$\delta=$' + str(delta) + r'$t=$' + str(t))

    ax.set_xlim(0,1)
    ax.set_xlabel(r'$x$')
    ax.set_ylabel(r'$P(x)$')
    ax.legend()
    
    if save_plot:
        folder = os.path.abspath('./plots/density/')
        if not os.path.exists(folder):
            os.makedirs(folder)
        file_name = os.path.join(folder, 'joint_ell{}_kappa{}_delta{}_t{}.pdf'.format(ell, kappa, delta, t))
        plt.savefig(file_name, bbox_inches='tight')
        plt.close()
    if show_plot:
        plt.show()



def read_csv(filename, model='v1'):
    """
    Read csv output for plotting
    """
    if model == 'v1':
        folder = os.path.abspath('./analytical_data/')
    elif model == 'v2':
        folder = os.path.abspath('./analytical_data_stop_start/')
    file = os.path.join(folder, filename)
    if not os.path.exists(file):
        raise ValueError("File does not exist: " + file)
    x_plot = []
    y_plot = []
    try:
        with open(file) as f:
            reader = csv.reader(f, delimiter=",")
            for line in reader:
                x_plot.append(float(line[0]))
                y_plot.append(float(line[1]))
    except:
        raise ValueError("Error reading file: " + file)
        
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

def plot_compare(ell, kappa, delta, t, phase, bins=20, ax=None):
    """
    Plot analytical and simulation data
    """
    filename = 'P_ell{}_kappa{}_delta{}_t{}.csv'.format(phase, ell, kappa, delta, t)
    x_plot, y_plot = read_csv(filename)
    if ax is None:
        fig, ax = initalize_plotting()
    ax.plot(x_plot, y_plot, label='Analytical')

    pos_file = get_pos_file(ell, kappa, delta, t, phase)
    x_all = read_pos(pos_file)
    prob = get_prob_phase(ell, kappa, delta, t, phase)
    
    hist, bin_edges = np.histogram(x_all, bins=bins, density=True)
    bin_centers = 0.5*(bin_edges[1:]+bin_edges[:-1])
    ax.plot(bin_centers, prob*hist, 'o', label='Simulation')
    ax.set_xlim(0,1)
    ax.set_xlabel(r'$x$')
    ax.set_ylabel(r'$P(x)$')
    ax.legend()
    
    return ax