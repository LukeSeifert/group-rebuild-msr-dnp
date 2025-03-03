import numpy as np
import matplotlib.pyplot as plt
import openmc.deplete
import openmc
import plotvals


def conc_collect(results, nucs):
    conc_dict = {}
    for nuc in nucs:
        conc_dict[nuc] = {}
        data = conc_dict[nuc]
        for r in results:
            data[r] = {}
            res = openmc.deplete.Results(r)
            t, conc = res.get_atoms('1', nuc, 'atom/cm3', 's')
            data[r]['x'] = t
            data[r]['y'] = conc
    return conc_dict

        
def plot_concs(conc_dict, names):
    markers = ['^', 'v', '<', '>']
    for nuc, data in conc_dict.items():
        res_iter = 0
        for key, value in data.items():
            plt.plot(value['x'], value['y'], label=names[res_iter], marker=markers[res_iter%len(markers)], markersize=5)
            res_iter += 1
        plt.xlabel(r'Time $[s]$')
        plt.ylabel(r'Concentration $[atoms/cm^3]$')
        plt.legend()
        plt.tight_layout()
        plt.savefig(f'concs_{nuc}.png')
        plt.close()
    return


if __name__ == '__main__':
    nucs = ['Br87']
    res_files = ['5.0-5.0times', '10.0-10.0times', '15.0-15.0times', '20.0-20.0times']
    names = [r'$\tau_{in}=5s, \tau_{ex}=5s$', r'$\tau_{in}=10s, \tau_{ex}=10s$',
             r'$\tau_{in}=15s, \tau_{ex}=15s$', r'$\tau_{in}=20s, \tau_{ex}=20s$'
             ]

    results = []
    concs_csv = []
    for r in res_files:
        res_path = f'./archived-data/results-{r}/Static/depletion_results.h5'
        conc_path = f'./archived-data/results-{r}/Static/concs.csv'
        results.append(res_path)
        concs_csv.append(conc_path)
        conc_dict = conc_collect(results, nucs)
    plot_concs(conc_dict, names)

    print('Done')