import matplotlib.pyplot as plt

nuc_list = ['Xe135']

default_batches = 100
default_nps = 1000
default_source = 1e10
default_temperature = 298.15
default_energy = 0.0253 * 1e-6
default_photons = False
default_run_mode = 'fixed source'
default_fissile = 'U235'
default_dens = 10
default_chain = '../../data/chain/chain_casl_pwr.xml'
default_incore_t = 4
default_excore_t = 4
default_name = 'default'
default_path = f'./results/{default_name}'


def dict_builder(name=default_name, incore_t=default_incore_t,
                 excore_t=default_excore_t, n_MeV=default_energy,
                 S_rate_per_s=default_source, batches=default_batches,
                 nps=default_nps, photon_bool=default_photons, 
                 run_mode=default_run_mode, temperature_K=default_temperature,
                 fissile_nuc=default_fissile, dens_g_cc=default_dens,
                 chain=default_chain):
    data = {
    'name': name,
    't_incore_s': incore_t,
    't_excore_s': excore_t,
    'n_MeV': n_MeV,
    'S_rate_per_s': S_rate_per_s,
    'batches': batches,
    'output_path': f'./results/{name}',
    'nps': nps,
    'photon_bool': photon_bool,
    'run_mode': run_mode,
    'temperature_K': temperature_K,
    'fissile_nuc': fissile_nuc,
    'dens_g_cc': dens_g_cc,
    'chain': chain
    }
    return data

base_case_data = dict_builder(name='flowing', incore_t=10, excore_t=10)
static_data = dict_builder(name='static', incore_t=20, excore_t=0,
                           S_rate_per_s=default_source/2)



















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
