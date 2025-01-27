from run import combine_pulse, run_fit
import itertools
import nlls
from copy import deepcopy
import shutil

def combinations(d):
    keys, values = zip(*d.items())
    return [dict(zip(keys, v)) for v in itertools.product(*values)]

if __name__ == '__main__':
    # Set up list of variables to change
    # Loop, move files, and run
    import ui

    run_omc = True
    decay_daughter = True

    all_fits = dict()
    all_fits['yield'] = {}
    all_fits['halflife'] = {}

    dt = 0.1
    tf = ui.default_omc_decay_time

    change_variables = {
        'nps': [10, 100, 500, 1000]
        #'temperature_K': [0, 250, 294, 600, 900, 1200, 2500]
    }

    new_vars = combinations(change_variables)
    pulse_data = deepcopy(ui.pulse_data)
    static_data = deepcopy(ui.static_data)

    for var_combo in new_vars:
        csv_name = str(list(var_combo.values())).strip('[]').strip(']').replace(', ', '-')

        pulse_data.update(var_combo)
        static_data.update(var_combo)

        all_fits = run_fit(all_fits, dt, tf, run_omc, decay_daughter,
                        ui.pulse_data, 'pulse')

        all_fits = run_fit(all_fits, dt, tf, run_omc, decay_daughter,
                        ui.static_data, 'saturation')

        all_fits = combine_pulse(all_fits)

        nlls.generate_csvs(all_fits, csv_name=csv_name)

        source = './results'
        destination = f'./postprocess/archived-data/results-{csv_name}'

        shutil.move(source, destination)