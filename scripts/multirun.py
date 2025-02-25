from run import combine_pulse, run_fit
import itertools
import nlls
from copy import deepcopy
import shutil
from time import time
import numpy as np

def combinations(d):
    keys, values = zip(*d.items())
    return [dict(zip(keys, v)) for v in itertools.product(*values)]

if __name__ == '__main__':
    # Set up list of variables to change
    # Loop, move files, and run
    import ui

    rate = 1/20
    repr_no_long = {'Kr': rate,
                'Xe': rate,
                'Se': rate,
                'Nb': rate,
                'Mo': rate,
                'Tc': rate,
                'Ru': rate,
                'Rh': rate,
                'Pd': rate,
                'Ag': rate,
                'Sb': rate,
                'Te': rate,
                }

    repr_long = {'Kr': rate,
                'Xe': rate,
                'Se': rate,
                'Nb': rate,
                'Mo': rate,
                'Tc': rate,
                'Ru': rate,
                'Rh': rate,
                'Pd': rate,
                'Ag': rate,
                'Sb': rate,
                'Te': rate,
                }

    rate = 1/(50*24*3600)
    more_data = {'Y': rate,
                'La': rate,
                'Ce': rate,
                'Pr': rate,
                'Nd': rate,
                'Pm': rate,
                'Sm': rate,
                'Gd': rate,
                'Eu': rate
    }
    repr_long.update(more_data)

    rate = 1/(60*24*3600)
    more_data = {'Br': rate,
                'I': rate
    }
    repr_long.update(more_data)

    rate = 1/(200*24*3600)
    more_data = {'Zr': rate,
                'Cd': rate,
                'In': rate,
                'Sn': rate
    }
    repr_long.update(more_data)



    run_omc = True
    decay_daughter = True
    num_times = 20
    repr_mults = [0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5]

    dt = 0.1
    tf = ui.default_omc_decay_time

    base_repr = repr_long.copy()
    multi_dict_eval = list()
    for mult in repr_mults:
        new_dict = {key: value * mult for key, value in base_repr.items()}
        multi_dict_eval.append(new_dict)

    change_variables = {
        #'nps': [1, 10, 100, 500, 1000, 5000, 10000, 50000, 100000, 500000, 1000000]
        #'temperature_K': [250, 294, 600, 900, 1200, 2500]
        #'final_time': [60, 120, 240, 420, 600] # time sample is irradiated
        #'dens_g_cc': [1, 5, 10, 50, 100]
        #'omc_dec_step': [0.1, 0.5, 1, 2, 5, 10]
        #'omc_dec_time': [60, 120, 240, 420, 600] # time sample is measured
        #'repr': [repr_no_long, repr_long]
        #
        'repr': multi_dict_eval # list of dicts with varying scaling rates
        #'t_incore_s': np.linspace(1, 20, num_times)
        #'t_excore_s': np.linspace(1, 20, num_times)
    }

    nameset = {
        'nps': 'nps',
        'temperature_K': 'K',
        'final_time': 'tirrad',
        'dens_g_cc': 'gpcc',
        'omc_dec_step': 'decdt',
        'omc_dec_time': 'dect',
        'repr': 'repr',
        't_incore_s': 'tin',
        't_excore_s': 'tex'
    }

    dopulse = {
        'nps': True,
        'temperature_K': True,
        'final_time': False,
        'dens_g_cc': True,
        'omc_dec_step': True,
        'omc_dec_time': True,
        'repr': True,
        't_incore_s': False,
        't_excore_s': False
    }

    if len(change_variables.keys()) == 1:
        naming_modifier = nameset[list(change_variables.keys())[0]]
        change_pulse_too = dopulse[list(change_variables.keys())[0]]
    elif len(change_variables.keys()) == 0:
        naming_modifier = 'daughter'
        change_pulse_too = True
        decay_daughter = False
    else:
        naming_modifier = 'times'
        change_pulse_too = False

    try:
        new_vars = combinations(change_variables)
    except ValueError:
        new_vars = [{}]
    pulse_data = deepcopy(ui.pulse_data)
    static_data = deepcopy(ui.static_data)
    time_taken = list()


    input(f'Change pulse is set to: {change_pulse_too} and name is {naming_modifier}')

    for combo_i, var_combo in enumerate(new_vars):
        start = time()
        all_fits = dict()
        all_fits['yield'] = {}
        all_fits['halflife'] = {}

        csv_name = str(list(var_combo.values())).strip('[]').strip(']').replace(', ', '-').replace('np.float64', '').replace('(', '').replace(')', '') + naming_modifier
        if csv_name == "{'Kr': 0.05-'Xe': 0.05-'Se': 0.05-'Nb': 0.05-'Mo': 0.05-'Tc': 0.05-'Ru': 0.05-'Rh': 0.05-'Pd': 0.05-'Ag': 0.05-'Sb': 0.05-'Te': 0.05}repr":
            csv_name = 'nolong'
        elif csv_name == "{'Kr': 0.05-'Xe': 0.05-'Se': 0.05-'Nb': 0.05-'Mo': 0.05-'Tc': 0.05-'Ru': 0.05-'Rh': 0.05-'Pd': 0.05-'Ag': 0.05-'Sb': 0.05-'Te': 0.05-'Y': 2.3148148148148148e-07-'La': 2.3148148148148148e-07-'Ce': 2.3148148148148148e-07-'Pr': 2.3148148148148148e-07-'Nd': 2.3148148148148148e-07-'Pm': 2.3148148148148148e-07-'Sm': 2.3148148148148148e-07-'Gd': 2.3148148148148148e-07-'Eu': 2.3148148148148148e-07-'Br': 1.9290123456790122e-07-'I': 1.9290123456790122e-07-'Zr': 5.787037037037037e-08-'Cd': 5.787037037037037e-08-'In': 5.787037037037037e-08-'Sn': 5.787037037037037e-08}repr":
            csv_name = 'long'
        elif len(csv_name) > 20 and naming_modifier == 'repr':
            csv_name = str(repr_mults[combo_i]) + 'repr'

        if change_pulse_too:
            pulse_data.update(var_combo)
        static_data.update(var_combo)

        all_fits = run_fit(all_fits, dt, tf, run_omc, decay_daughter,
                        pulse_data, 'pulse')

        all_fits = run_fit(all_fits, dt, tf, run_omc, decay_daughter,
                        static_data, 'saturation')

        all_fits = combine_pulse(all_fits)

        nlls.generate_csvs(all_fits, csv_name=csv_name)

        source = './results'
        destination = f'./postprocess/archived-data/results-{csv_name}'

        shutil.move(source, destination)
        end = time()
        time_taken.append(end-start)
    
    for ti, t in enumerate(time_taken):
        print(f'Number {ti+1} took {round(t, 2)} s')
    print(f'Net time of {round(sum(time_taken), 3)} s')