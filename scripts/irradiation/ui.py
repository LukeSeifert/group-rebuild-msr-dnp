import matplotlib.pyplot as plt

name = 'flowing'
base_case_data = {
'name': name,
't_incore_s': 10,
't_excore_s': 10,
'n_MeV': 0.0253 * 1e-6,
'S_rate_per_s': 1e10,
'batches': 10,
'output_path': f'./results/{name}',
'nps': 1000,
'photon_bool': False,
'run_mode': 'fixed source',
'temperature_K': 298.15,
'fissile_nuc': 'U235',
'dens_g_cc': 10,
'chain': '../../data/chain/chain_casl_pwr.xml'
}

name = 'static'
static_data = {
'name': name,
't_incore_s': 20,
't_excore_s': 0,
'n_MeV': 0.0253 * 1e-6,
'S_rate_per_s': 1e10,
'batches': 10,
'output_path': f'./results/{name}',
'nps': 1000,
'photon_bool': False,
'run_mode': 'fixed source',
'temperature_K': 298.15,
'fissile_nuc': 'U235',
'dens_g_cc': 10,
'chain': '../../data/chain/chain_casl_pwr.xml'
}




















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
