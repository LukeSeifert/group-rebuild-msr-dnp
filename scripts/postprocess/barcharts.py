import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams["text.usetex"] = "True"
plt.rcParams["font.size"] = 16
plt.rcParams["axes.labelsize"] = 20
plt.rcParams["axes.labelweight"] = "bold"
plt.rcParams["lines.linewidth"] = 1.5
plt.rcParams["lines.markersize"] = 1
plt.rcParams["axes.grid"] = True
plt.rcParams["axes.grid.which"] = "major"
plt.rcParams["grid.linestyle"] = "--"
plt.rcParams["grid.linewidth"] = 1
plt.rcParams["xtick.direction"] = "in"
plt.rcParams["ytick.direction"] = "in"
plt.rcParams["xtick.major.size"] = 6.0
plt.rcParams["ytick.major.size"] = 6.0
plt.rcParams["xtick.minor.size"] = 3.0
plt.rcParams["ytick.minor.size"] = 3.0
plt.rcParams["figure.autolayout"] = True
plt.rcParams['savefig.dpi'] = 300

data_name = '560cp3'
plot_topics = ['yields', 'halflives']
y_names = ['Yield', 'Half-life [s]']
yields = dict()
halflives = dict()
for i, fname in enumerate(plot_topics):
    try:
        df = pd.read_csv(f'./{fname}/{data_name}.csv')
    except FileNotFoundError:
        print(f'{fname} not available')
        continue

    if fname == 'yields':
        data_sources = list(set(df["Data Source"].tolist()))
        for data_source in data_sources:
            yields[data_source] = df.loc[df["Data Source"] == data_source, f"{y_names[i]}"]
            net_yield = yields[data_source].sum()
            print(f'{data_source} yield: {round(net_yield, 5)}')
    
    if fname == 'halflives':
        plt.yscale('log')
        data_sources = list(set(df["Data Source"].tolist()))
        for data_source in data_sources:
            halflives[data_source] = df.loc[df["Data Source"] == data_source, f"{y_names[i]}"]
            avg_hl = np.sum(yields[data_source] * halflives[data_source] / yields[data_source].sum())
            print(f'{data_source} average half-life: {round(avg_hl, 3)} s')

    sns.barplot(df, x='Group', y=y_names[i], hue='Data Source')
    plt.legend(fontsize=12)
    plt.savefig(f'{fname}.png')
    plt.close()
