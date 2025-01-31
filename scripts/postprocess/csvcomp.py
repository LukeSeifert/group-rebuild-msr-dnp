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
    csv_name = './post-data.csv'
    yields = dict()
    halflives = dict()
    names = list()
    net_yields = list()
    avg_hls = list()
    for i, fname in enumerate(plot_topics):
        df_use = None
        for j, csv in enumerate(csvs):
            #print(f'\nFile: {csv}')
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
                df_new[compare_column] = df_new[compare_column].replace(compare_column_value, csv_name[j])
            try:
                df_use = pd.concat([df_use, df_new])
            except TypeError:
                df_use = df_new

        try:
            sns.barplot(df_use, x='Group', y=y_names[i], hue='Data Source')
            plt.legend(fontsize=12)
            plt.savefig(f'{fname}-csvs.png')
            plt.close()
        except ValueError:
            continue
    data = dict()
    data['Data Source'] = names
    data['Yield'] = net_yields
    data['Half-life [s]'] = avg_hls
    df = df.from_dict(data)
    print(df)
    df.to_csv('yield-hl-data.csv')
    return

def static_concentration_comparison(csvs, top_num=10):
    running_df = None
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
    for nuc in all_nucs:
        cur_nuc_rows = running_df.loc[running_df['Nuclide'] == nuc]
        cur_nuc_concs = cur_nuc_rows['Concentration']
        std_dev = np.std(cur_nuc_concs)
        mean = np.mean(cur_nuc_concs)
        std_dev_concs[nuc] = std_dev / mean # coefficient of variation
    
    sorted_devs = sorted(zip(std_dev_concs.values(), std_dev_concs.keys()),
                         reverse=True)[:top_num]
    for std_dev, nuc in sorted_devs:
        print(f'{nuc} - {round(std_dev*100, 3)}%')
        sns.barplot(running_df, x='Parameter', y='Concentration', errorbar=None)
        plt.savefig(f'csvs-{nuc}.png')
        plt.close()
    return
        
        



if __name__ == '__main__':
    nps_analysis = True
    temp_analysis = False
    num_nucs = 1

    if nps_analysis:
        csvs = ['1nps', '10nps', '100nps', '500nps', '1000nps', '5000nps', '10000nps']
        csv_name = ['1', '10', '100', '500', '1000', '5000', '10000']
    elif temp_analysis:
        csvs = ['0K', '250K', '294K', '600K', '900K', '920K', '1200K', '2500K']
        csv_name = ['0K', '250K', '294K', '600K', '900K', '920K', '1200K', '2500K']
    else:
        raise Exception('No analysis selected')

    static_concentration_comparison(csvs, top_num=num_nucs)
    plot_comparison(csvs, csv_name)