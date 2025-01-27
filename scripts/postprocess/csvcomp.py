import matplotlib.pyplot as plt
import plotvals
from barcharts import collect_data
import seaborn as sns
import pandas as pd
import warnings

csvs = ['default', '920K']
csv_name = ['298.15K', '920K']
compare_column = 'Data Source'
compare_column_value = 'Static-Pulse'
plot_topics = ['yields', 'halflives']
y_names = ['Yield', 'Half-life [s]']
yields = dict()
halflives = dict()
for i, fname in enumerate(plot_topics):
    df_use = None
    for j, csv in enumerate(csvs):
        print(f'\nFile: {csv}')
        df, yields, halflives, net_yield, avg_hl = collect_data(fname,
                                                                csv,
                                                                y_names,
                                                                yields,
                                                                halflives,
                                                                i)
        df_new = df[df[compare_column] == compare_column_value]
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
