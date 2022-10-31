from matplotlib import pyplot as plt
import json
import numpy as np
from cycler import cycler
import matplotlib.colors as mcolors

# Load the dictionary
def load_data(filename):
    with open(filename) as json_file:
        return json.load(json_file)

# Get the pareto data
def pareto_data(filename):
    dict = load_data(filename)
    alphas = []
    costs = []
    for k, v in dict.items():
        if v['status'] == 'optimal':
            alphas.append(float(k))
            costs.append(float(v['objective']))
    return costs, alphas

# Plot the Pareto tradeoff
def plot_pareto_comparison(files, save=False):
    data = {}
    max_cost = 0
    for label, file in files.items():
        costs, alphas = pareto_data(file)
        data[label] = (costs, alphas)
        highest = max(costs)
        if highest > max_cost:
            max_cost = highest
    for v in data.values():
        v[0].append(max_cost)
        v[1].append(v[1][-1])

    fig, ax = plt.subplots(1,1)
    custom_cycler = (cycler(marker=["o", "s", "^", "p", "D", "*", "8"]) +
                     cycler(color=mcolors.to_rgba_array(['C0', 'C1', 'C2', 'C3', 'C4', 'C5', 'C6'])))
    ax.set_prop_cycle(custom_cycler)
    plt.xlabel('Cost')
    plt.ylabel(r'Probability Level $\alpha$')
    for label, (costs, alphas) in data.items():
        plt.plot(costs, alphas, label=label)
    plt.legend(loc='best')
    plt.ylim(0.8,1.01)
    if save:
        plt.savefig("pareto.png", dpi=300, transparent=True)
    plt.show()

# Get the time data
def time_data(file):
    dict = load_data(file)
    return np.sort(np.array([float(v['time']) for v in dict.values()]))

# Plot the performance data 
def plot_performance(files, save=False):
    data = {}
    max_time = 0.0
    for label, file in files.items():
        ts = time_data(file)
        if ts[-1] > max_time:
            max_time = ts[-1]
        fractions = [(i + 1) / len(ts) for i in range(len(ts)) if ts[i] < 3600]
        if len(fractions) < len(ts):
            fractions += (len(ts) - len(fractions)) * [fractions[-1]]
        data[label] = (ts.tolist(), fractions)
    for v in data.values():
        v[0].append(max_time)
        v[1].append(v[1][-1])

    fig, ax = plt.subplots(1,1)
    custom_cycler = (cycler(linestyle=["-", "--", ":", "-.", (0, (3, 5, 1, 5, 1, 5))]) +
                     cycler(color=mcolors.to_rgba_array(['C0', 'C1', 'C2', 'C3', 'C4'])))
    ax.set_prop_cycle(custom_cycler)
    plt.xlabel('Time (s)')
    plt.ylabel('Fraction of Instances Solved')
    for label, v in data.items():
        plt.plot(v[0], v[1], label=label)
    plt.legend(loc='best')
    plt.xlim(0, 3700)
    if save:
        plt.savefig("performance.png", dpi=300, transparent=True)
    plt.show()

# Time Comparison
def plot_time_comparison(files, save=False):
    xs = []
    plt.figure()
    fig, ax = plt.subplots()
    x = 0
    labels = [20,19,18,17,16,15]*3

    for label, file in files.items():
        xs.append(x)
        ts = time_data(file)
        avg = np.mean(ts)
        std = np.std(ts)/2
        color = None
        hatch = None
        if label[0] == '5':
            color = 'C0'
            hatch = '///'
        elif label[0] == '4':
            color = 'C1'
            hatch = '\\\\\\'
        elif label[0] == '3':
            color = 'C2'
            hatch = 'xxx'
        
        if label == '5 gens and 19 lines':
            ax.bar(x, avg, yerr=std, color=color, capsize=5, hatch=hatch, alpha = 0.6, label=r'$Y_{g, min} = 5$')
        elif label == '4 gens and 19 lines':
            ax.bar(x, avg, yerr=std,color=color, capsize=5, hatch=hatch, alpha = 0.6, label=r'$Y_{g, min} = 4$')
        elif label == '3 gens and 19 lines':
            ax.bar(x, avg, yerr=std, color=color, capsize=5, hatch=hatch, alpha = 0.6, label=r'$Y_{g, min} = 3$')
        else:
            ax.bar(x, avg, yerr=std, capsize=5, hatch=hatch, alpha = 0.6, color=color)
        x += 1
        
    plt.xticks(xs, labels)
    plt.legend(loc='best')
    plt.xlabel(r'$Y_{l, min}$')
    plt.ylabel('Average Time (s)')
    if save: 
        plt.savefig("time.png", dpi=300, transparent=True)
    plt.show()

# Make the plots
labeled_files_pareto = {r'$\Omega_\wedge(Y_g, Y_l)$' : './data/ieee14_intersection_bigm_1000.json',
                   r'$\Omega_{atleast}(Y_g, Y_l; 5, 19)$' : './data/ieee14_atleast_5gens_19lines_bigm_1000.json',
                   r'$\Omega_{atleast}(Y_g, Y_l; 5, 18)$' : './data/ieee14_atleast_5gens_18lines_bigm_1000.json',
                   r'$\Omega_{atleast}(Y_g, Y_l; 4, 20)$' : './data/ieee14_atleast_4gens_20lines_bigm_1000.json',
                   r'$\Omega_{atleast}(Y_g, Y_l; 3, 15)$' : './data/ieee14_atleast_3gens_15lines_bigm_1000.json',
                   }

labeled_files_time = {'5 gens and 20 lines' : './data/ieee14_intersection_bigm_1000.json',
                   '5 gens and 19 lines' : './data/ieee14_atleast_5gens_19lines_bigm_1000.json',
                   '5 gens and 18 lines' : './data/ieee14_atleast_5gens_18lines_bigm_1000.json',
                   '5 gens and 17 lines' : './data/ieee14_atleast_5gens_17lines_bigm_1000.json',
                   '5 gens and 16 lines' : './data/ieee14_atleast_5gens_16lines_bigm_1000.json',
                   '5 gens and 15 lines' : './data/ieee14_atleast_5gens_15lines_bigm_1000.json',
                   '4 gens and 20 lines' : './data/ieee14_atleast_4gens_20lines_bigm_1000.json',
                   '4 gens and 19 lines' : './data/ieee14_atleast_4gens_19lines_bigm_1000.json',
                   '4 gens and 18 lines' : './data/ieee14_atleast_4gens_18lines_bigm_1000.json',
                   '4 gens and 17 lines' : './data/ieee14_atleast_4gens_17lines_bigm_1000.json',
                   '4 gens and 16 lines' : './data/ieee14_atleast_4gens_16lines_bigm_1000.json',
                   '4 gens and 15 lines' : './data/ieee14_atleast_4gens_15lines_bigm_1000.json',
                   '3 gens and 20 lines' : './data/ieee14_atleast_3gens_20lines_bigm_1000.json',
                   '3 gens and 19 lines' : './data/ieee14_atleast_3gens_19lines_bigm_1000.json',
                   '3 gens and 18 lines' : './data/ieee14_atleast_3gens_18lines_bigm_1000.json',
                   '3 gens and 17 lines' : './data/ieee14_atleast_3gens_17lines_bigm_1000.json',
                   '3 gens and 16 lines' : './data/ieee14_atleast_3gens_16lines_bigm_1000.json',
                   '3 gens and 15 lines' : './data/ieee14_atleast_3gens_15lines_bigm_1000.json',
                   }


labeled_files_perform = {'Big-M' : './data/ieee14_intersection_bigm_1000.json',
                         'Hull' : './data/ieee14_intersection_hull_1000.json',
                         'Cutting-Plane' : './data/ieee14_intersection_cuttingplane_1000.json'}#,
                        #  'RIC' : './data/ieee14_intersection_RIC_1000.json',
                        #  'LOA' : './data/ieee14_intersection_LOA_1000.json'}

plot_pareto_comparison(labeled_files_pareto, True)

plot_performance(labeled_files_perform, True)

plot_time_comparison(labeled_files_time, True)
