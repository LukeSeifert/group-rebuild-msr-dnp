import matplotlib.pyplot as plt

nuc_list = ['Rb91', 'Br87', 'Cs141']

default_batches = 10
default_nps = 100
default_source = 4.2e16
default_temperature = 298.15
default_energy = 0.0253 * 1e-6
default_photons = False
default_run_mode = 'fixed source'
default_fissile = 'U235'
default_dens = 10
default_incore_t = 4
default_excore_t = 4
default_name = 'Default'
default_path = f'./results/{default_name}'
default_final_time = 7 * 60

if default_energy < 300 * 1e-6:
    default_chain = '../data/chain/chain_endfb80_pwr.xml'
else:
    default_chain = '../data/chain/chain_endfb80_sfr.xml'




def dict_builder(name=default_name, incore_t=default_incore_t,
                 excore_t=default_excore_t, n_MeV=default_energy,
                 net_Sp=default_source, batches=default_batches,
                 nps=default_nps, photon_bool=default_photons, 
                 run_mode=default_run_mode, temperature_K=default_temperature,
                 fissile_nuc=default_fissile, dens_g_cc=default_dens,
                 chain=default_chain, final_time=default_final_time):
    S_rate_per_s = net_Sp * (excore_t + incore_t) / (incore_t * final_time)
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
    'chain': chain,
    'final_time': final_time
    }
    return data

flow_data = dict_builder(name='Flowing', incore_t=10, excore_t=10)

static_data = dict_builder(name='Static', incore_t=10*2, excore_t=0)

pulse_data = dict_builder(name='Pulse', incore_t=1/4*1e-3, excore_t=0,
                          final_time=10/4*1e-3)

mostly_excore_data = dict_builder(name='ExFlowing', incore_t=10/10,
                                  excore_t=10)



















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
