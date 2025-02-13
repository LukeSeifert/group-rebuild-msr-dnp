import matplotlib.pyplot as plt
import plotvals
from analysis import collect_data
import seaborn as sns
import pandas as pd
import warnings
import numpy as np

def plot_comparison(csvs, csv_name):
    compare_column = 'Data Source'
    compare_column_value = 'Static-Pulse'
    plot_topics = ['yields', 'halflives']
    y_names = ['Yield', 'Half-life [s]']
    yields = dict()
    halflives = dict()
    names = list()
    net_yields = list()
    avg_hls = list()
    for i, fname in enumerate(plot_topics):
        df_use = None
        for j, csv in enumerate(csvs):
            df, yields, halflives, net_yield, avg_hl = collect_data(fname,
                                                                    csv,
                                                                    y_names,
                                                                    yields,
                                                                    halflives,
                                                                    i,
                                                                    compare_column_value=compare_column_value)
            if net_yield > 0:
                names.append(csv)
                net_yields.append(net_yield)
            if avg_hl > 0:
                avg_hls.append(avg_hl)

            try:
                df_new = df[df[compare_column] == compare_column_value]
            except TypeError:
                continue
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                df_new[compare_column] = df_new[compare_column].replace(compare_column_value, csv)
                df_new['Stripped Source'] = csv_name[j]
            try:
                df_use = pd.concat([df_use, df_new])
            except TypeError:
                df_use = df_new

        try:
            sns.barplot(df_use, x='Group', y=y_names[i], hue='Stripped Source')
            plt.legend(fontsize=12)
            plt.savefig(f'{fname}-csvs.png')
            plt.close()
        except ValueError:
            continue
    data = dict()
    data['Data Source'] = csv_name
    data[r'$\bar{\nu}_d$'] = net_yields
    data[r'$\bar{T}$'] = avg_hls
    data[r'$|\Delta \bar{\nu}_d|$'] = [abs(i - net_yields[-1]) for i in net_yields]
    data[r'$|\Delta \bar{T}| [s]$'] = [abs(i - avg_hls[-1]) for i in avg_hls]
    df = df.from_dict(data)
    print(df.to_latex(index=False))
    df.to_csv('yield-hl-data.csv')
    return

def static_concentration_comparison(csvs, top_num=10):
    running_df = None
    new_df_data = dict()
    new_df_data['Nuclide'] = list()
    new_df_data['CV [%]'] = list()
    new_df_data['Pn'] = list()
    new_df_data['Half-life [s]'] = list()
    new_df_data['Eval-term'] = list()
    for csv in csvs:
        df = pd.read_csv(f'./archived-data/results-{csv}/Static/concs.csv',
                         header=None)
        df.columns = ['Nuclide', 'Concentration']
        df['Parameter'] = csv
        try:
            running_df = pd.concat([running_df, df], ignore_index=True)
        except TypeError:
            running_df = df
    
    all_nucs = list(set(df['Nuclide']))
    std_dev_concs = dict()
    concentration_averages = dict()
    for nuc in all_nucs:
        cur_nuc_rows = running_df.loc[running_df['Nuclide'] == nuc]
        cur_nuc_concs = cur_nuc_rows['Concentration']
        std_dev = np.std(cur_nuc_concs)
        mean = np.mean(cur_nuc_concs)
        concentration_averages[nuc] = mean
        std_dev_concs[nuc] = std_dev / mean # coefficient of variation
    
    data_df = pd.read_csv(f'./archived-data/results-{csv}/Static/data.csv')
    net_dn_yield = 0
    for nuc in all_nucs:
        cur_data = data_df.loc[data_df['Nuclide'] == nuc]
        pn = float(cur_data['Pn'])
        lam = float(cur_data['lam'])
        dn_yield = pn*lam*concentration_averages[nuc]
        net_dn_yield += dn_yield

    for nuc in all_nucs:
        cur_data = data_df.loc[data_df['Nuclide'] == nuc]
        pn = float(cur_data['Pn'])
        lam = float(cur_data['lam'])
        new_df_data['Nuclide'].append(nuc)
        new_df_data['CV [%]'].append(std_dev_concs[nuc] * 100)
        new_df_data['Pn'].append(pn)
        new_df_data['Half-life [s]'].append(np.log(2)/lam)
        dn_yield = pn*lam*concentration_averages[nuc]
        new_df_data['Eval-term'].append(dn_yield/net_dn_yield*std_dev_concs[nuc])

    new_df = pd.DataFrame.from_dict(new_df_data)
    new_df = new_df.sort_values(by='Eval-term', ignore_index=True,
                                ascending=False).iloc[:top_num]
    print(new_df)
    return
        
        



if __name__ == '__main__':
    nps_analysis = False
    temp_analysis = False
    tirrad_analysis = False
    dens_analysis = False
    decdt_analysis = False
    dect_analysis = False
    longrepr_analysis = True
    num_nucs = 10

    if nps_analysis:
        csvs = ['1nps', '10nps', '100nps', '500nps', '1000nps', '5000nps', '10000nps', '50000nps', '100000nps', '500000nps', '1000000nps']
        csv_name = ['1', '10', '100', '500', '1000', '5000', '10000', '50000', '100000', '500000', '1000000']
    elif temp_analysis:
        csvs = ['250K', '294K', '600K', '900K', '920K', '1200K', '2500K']
        csv_name = ['250K', '294K', '600K', '900K', '920K', '1200K', '2500K']
    elif tirrad_analysis:
        csvs =  ['60tirrad', '120tirrad', '240tirrad', '420tirrad', '600tirrad']
        csv_name = [r'$60s$', r'$120s$', r'$240s$', r'$420s$', r'$600s$']
    elif dens_analysis:
        csvs = ['1gpcc', '5gpcc', '10gpcc', '50gpcc', '100gpcc']
        csv_name = [r'$1\frac{g}{cm^3}$', r'$5\frac{g}{cm^3}$',
                    r'$10\frac{g}{cm^3}$', r'$50\frac{g}{cm^3}$',
                    r'$100\frac{g}{cm^3}$']
    elif decdt_analysis:
        base = [0.1, 0.5, 1, 2, 5, 10]
        base.reverse()
        csvs = [str(i) + 'decdt' for i in base]
        csv_name = [str(i) + r'$s$' for i in base]
    elif dect_analysis:
        base = [60, 120, 240, 420, 600]
        csvs = [str(i) + 'dect' for i in base]
        csv_name = [str(i) + r'$s$' for i in base]
    elif longrepr_analysis:
        csvs = ['nolong', 'long']
        csv_name = ['Without Large', 'With Large']
    else:
        raise Exception('No analysis selected')

    static_concentration_comparison(csvs, top_num=num_nucs)
    plot_comparison(csvs, csv_name)