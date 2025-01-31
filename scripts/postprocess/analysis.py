import seaborn as sns
import plotvals
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


def collect_data(fname, data_name, y_names, yields, halflives, i,
                 data_sources:list=None, net_yield=0,
                 avg_hl=0):
    try:
        df = pd.read_csv(f'./{fname}/{data_name}.csv')
    except FileNotFoundError:
        print(f'{fname} not available')
        return [None] * 5

    if fname == 'yields':
        if type(data_sources) == type(None):
            data_sources = list(set(df["Data Source"].tolist()))
        for data_source in data_sources:
            yields[data_source] = df.loc[df["Data Source"] == data_source, f"{y_names[i]}"]
            net_yield = yields[data_source].sum()
            print(f'{data_source} yield: {round(net_yield, 5)}')
    
    if fname == 'halflives':
        plt.yscale('log')
        if type(data_sources) == type(None):
            data_sources = list(set(df["Data Source"].tolist()))
        for data_source in data_sources:
            halflives[data_source] = df.loc[df["Data Source"] == data_source, f"{y_names[i]}"]
            avg_hl = np.sum(yields[data_source] * halflives[data_source] / yields[data_source].sum())
            print(f'{data_source} average half-life: {round(avg_hl, 3)} s')
    
    return df, yields, halflives, net_yield, avg_hl


if __name__ == '__main__':
    data_name = '920K'
    plot_topics = ['yields', 'halflives']
    y_names = ['Yield', 'Half-life [s]']
    yields = dict()
    halflives = dict()
    for i, fname in enumerate(plot_topics):
        df, yields, halflives, net_yield, avg_hl = collect_data(fname,
                                                                data_name,
                                                                y_names,
                                                                yields,
                                                                halflives,
                                                                i)

        try:
            sns.barplot(df, x='Group', y=y_names[i], hue='Data Source')
            plt.legend(fontsize=12)
            plt.savefig(f'{fname}.png')
            plt.close()
        except ValueError:
            continue
