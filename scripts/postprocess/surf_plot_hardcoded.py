import plotvals
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from matplotlib import cm
from matplotlib.ticker import LinearLocator



plt.rcParams["grid.linestyle"] = ""

# Make data.
X = [5,
5,
5,
5,
10,
10,
10,
10,
15,
15,
15,
15,
20,
20,
20,
20]
Y = [5,
10,
15,
20,
5,
10,
15,
20,
5,
10,
15,
20,
5,
10,
15,
20]
Z = [
0.01594,
0.015165,
0.014179,
0.013929,
0.015851,
0.01479,
0.014109,
0.013845,
0.01622,
0.014718,
0.014051,
0.013693,
0.019775,
0.014672,
0.013996,
0.022849
]

mindex = Z.index(min(Z))
maxdex = Z.index(max(Z))
print(f'Min of {Z[mindex]} for in {X[mindex]} and ex {Y[mindex]}')
print(f'Max of {Z[maxdex]} for in {X[maxdex]} and ex {Y[maxdex]}')
print(f'Diff: {(Z[maxdex] - Z[mindex]) * 1e5}pcm')

x = np.reshape(X, (4, 4))
y = np.reshape(Y, (4, 4))
z = np.reshape(Z, (4, 4))

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
#ax = sns.heatmap(pivotted,cmap='cubehelix')
ax = sns.heatmap(pivotted, cmap=color)
ax.collections[0].colorbar.set_label(z_name)
plt.tight_layout()

#fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
#surf = ax.plot_surface(x, y, z, cmap=cm.coolwarm,
#                       linewidth=0, antialiased=False)

#ax.zaxis.set_major_locator(LinearLocator(10))
#ax.zaxis.set_major_formatter('{x:.02f}')
#ax.set_xlabel('In-core Residence Time [s]')
#ax.set_ylabel('Ex-core Residence Time [s]')
#ax.set_zlabel('Delayed Neutron Yield')
#fig.colorbar(surf, shrink=0.5, aspect=5)

plt.savefig('surf.png')
#plt.show()
plt.close()


plt.rcParams["grid.linestyle"] = "--"