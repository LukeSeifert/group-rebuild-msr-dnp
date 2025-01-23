import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd

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

data_name = 'fastu235ORIGEN'
plot_topics = ['yields', 'halflives']
y_names = ['Yield', 'Half-life [s]']
for i, fname in enumerate(plot_topics):
    try:
        df = pd.read_csv(f'./{fname}/{data_name}.csv')
    except FileNotFoundError:
        print(f'{fname} not available')
        continue

    if fname == 'yields':
        data_sources = list(set(df["Data Source"].tolist()))
        for data_source in data_sources:
            net_yield = df.loc[df["Data Source"] == data_source, f"{y_names[i]}"].sum()
            print(f'{data_source} yield: {round(net_yield, 5)}')
    
    if fname == 'halflives':
        plt.yscale('log')

    sns.barplot(df, x='Group', y=y_names[i], hue='Data Source')
    plt.legend(fontsize=12)
    plt.savefig(f'{fname}.png')
    plt.close()
