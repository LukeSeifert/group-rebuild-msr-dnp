import numpy as np
import matplotlib.pyplot as plt
import openmc.deplete
import openmc
import plotvals
from collections import OrderedDict


def conc_collect(results, nucs):
    conc_dict = OrderedDict()
    for ri, r in enumerate(results):
        for nuc in nucs:
            if ri == 0:
                conc_dict[nuc] = {}
            conc_dict[nuc][r] = {}
            res = openmc.deplete.Results(r)
            t, conc = res.get_atoms('1', nuc, 'atom/cm3', 's')
            conc_dict[nuc][r]['x'] = t
            conc_dict[nuc][r]['y'] = conc
    return conc_dict


def combine_dict(conc_dict, deca_dict):
    net_data = OrderedDict()
    for nuc in conc_dict.keys():
        net_data[nuc] = dict()
        for i, key in enumerate(conc_dict[nuc].keys()):
            net_data[nuc][key] = dict()
            conc_data = list(conc_dict[nuc].values())[i]
            deca_data = list(deca_dict[nuc].values())[i]
            net_data[nuc][key]['x'] = np.append(conc_data['x'], deca_data['x'] + conc_data['x'][-1])
            net_data[nuc][key]['y'] = np.append(conc_data['y'], deca_data['y'])
    return net_data

 
def plot_concs(conc_dict, names, typing='concs'):
    markers = ['^', 'v', '<', '>']
    for nuc, data in conc_dict.items():
        res_iter = 0
        for key, value in data.items():
            if typing == 'concs':
                plt.plot(value['x'], value['y'], label=names[res_iter], marker=markers[res_iter%len(markers)], markersize=4, linestyle='--',
                        markevery=0.1)
                plt.yscale('linear')
            else:
                plt.plot(value['x'], value['y'], label=names[res_iter], marker=markers[res_iter%len(markers)], markersize=4, linestyle='--',
                        markevery=0.1)
                plt.yscale('log')
            res_iter += 1
        plt.xlabel(r'Time $[s]$')
        plt.ylabel(r'Concentration $[atoms/cm^3]$')
        plt.legend()
        plt.tight_layout()
        plt.savefig(f'{typing}_{nuc}.png')
        plt.close()
    return


if __name__ == '__main__':
    nucs = ['Br87']#['Br87', 'As86', 'Ge86', 'Br90', 'I137', 'Rb95', 'Br91']
    res_files = ['5.0-0times', '25.0-0times']
    names = [r'$\tau_{in}=5s, \tau_{ex}=0s$', r'$\tau_{in}=25s, \tau_{ex}=0s$'
             ]
    #res_files = ['5.0-5.0times', '10.0-10.0times', '15.0-15.0times', '20.0-20.0times']
    #names = [r'$\tau_{in}=5s, \tau_{ex}=5s$', r'$\tau_{in}=10s, \tau_{ex}=10s$',
    #         r'$\tau_{in}=15s, \tau_{ex}=15s$', r'$\tau_{in}=20s, \tau_{ex}=20s$'
    #         ]

    results = []
    decays = []
    for r in res_files:
        res_path = f'./archived-data/results-{r}/Static/depletion_results.h5'
        decay_path = f'./archived-data/results-{r}/StaticDecay/depletion_results.h5'
        results.append(res_path)
        decays.append(decay_path)
    conc_dict = conc_collect(results, nucs)
    deca_dict = conc_collect(decays, nucs)
    net_dict = combine_dict(conc_dict, deca_dict)

    plot_concs(conc_dict, names, typing='concs')
    plot_concs(deca_dict, names, typing='decay')
    plot_concs(net_dict, names, typing='net')

    print('Done')