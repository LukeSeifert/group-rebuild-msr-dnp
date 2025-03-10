import plotvals
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from matplotlib import cm
from matplotlib.ticker import LinearLocator



plt.rcParams["grid.linestyle"] = ""

# Make data.
X = [
5, 5, 5, 5, 5, 5, 
10, 10, 10, 10, 10, 10,
15, 15, 15, 15, 15, 15,
20, 20, 20, 20, 20, 20,
25, 25, 25, 25, 25, 25,
30, 30, 30, 30, 30, 30
]
Y = [
0, 5, 10, 15, 20, 25,
0, 5, 10, 15, 20, 25,
0, 5, 10, 15, 20, 25,
0, 5, 10, 15, 20, 25,
0, 5, 10, 15, 20, 25,
0, 5, 10, 15, 20, 25
]
Z = [
0.018736,
0.017171,
0.018734,
0.020703,
0.022045,
0.024553,
0.018741,
0.017389,
0.018398,
0.019002,
0.02118,
0.022541,
0.018731,
0.017515,
0.017354,
0.0194,
0.020506,
0.018333,
0.018739,
0.016646,
0.018222,
0.019181,
0.01619,
0.014798,
0.018636,
0.017614,
0.024013,
0.014761,
0.013496,
0.014746,
0.018745,
0.017644,
0.014115,
0.013111,
0.013785,
0.017555
]

dims = 6
mindex = Z.index(min(Z))
maxdex = Z.index(max(Z))
print(f'Min of {Z[mindex]} for in {X[mindex]} and ex {Y[mindex]}')
print(f'Max of {Z[maxdex]} for in {X[maxdex]} and ex {Y[maxdex]}')
print(f'Diff: {(Z[maxdex] - Z[mindex]) * 1e5}pcm')

x = np.reshape(X, (dims, dims))
y = np.reshape(Y, (dims, dims))
z = np.reshape(Z, (dims, dims))

x_name = r'$\tau_{in}$ $[s]$ '
y_name = r'$\tau_{ex}$ $[s]$'
z_name = r'$\bar{\nu}_d$'

#x_name = r'In-core Residence Time $[s]$'
#y_name = r'Ex-core Residence Time $[s]$'
#z_name = 'Total Delayed Neutron Yield'


df = pd.DataFrame.from_dict(np.array([X,Y,Z]).T)
df.columns = [x_name,y_name,z_name]
df[z_name] = pd.to_numeric(df[z_name])
pivotted= df.pivot(columns=y_name,index=x_name,values=z_name)
color = sns.color_palette("dark:pink_r", as_cmap=True)
ax = sns.heatmap(pivotted, cmap=color)
ax.invert_yaxis()
ax.collections[0].colorbar.set_label(z_name)
plt.tight_layout()


plt.savefig('surf.png')
#plt.show()
plt.close()

#fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
#surf = ax.plot_surface(x, y, z, cmap=cm.coolwarm,
#                       linewidth=0, antialiased=False)
#
#ax.zaxis.set_major_locator(LinearLocator(10))
#ax.zaxis.set_major_formatter('{x:.02f}')
#ax.set_xlabel('In-core Residence Time [s]')
#ax.set_ylabel('Ex-core Residence Time [s]')
#ax.set_zlabel('Delayed Neutron Yield')
#fig.colorbar(surf, shrink=0.5, aspect=5)
#plt.show()
#plt.close()

plt.rcParams["grid.linestyle"] = "--"