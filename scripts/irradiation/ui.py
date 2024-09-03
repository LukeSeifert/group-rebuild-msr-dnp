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