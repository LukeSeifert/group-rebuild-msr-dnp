from run import combine_pulse, run_fit
import itertools
import nlls
from copy import deepcopy
import shutil
from time import time

def combinations(d):
    keys, values = zip(*d.items())
    return [dict(zip(keys, v)) for v in itertools.product(*values)]

if __name__ == '__main__':
    # Set up list of variables to change
    # Loop, move files, and run
    import ui

    naming_modifier = 'nps'

    run_omc = True
    decay_daughter = True

    dt = 0.1
    tf = ui.default_omc_decay_time

    change_variables = {
        'nps': [5000, 10000]
        #'nps': [1, 10, 100, 500, 1000]
        #'temperature_K': [250, 294, 600, 900, 1200, 2500]
    }

    new_vars = combinations(change_variables)
    pulse_data = deepcopy(ui.pulse_data)
    static_data = deepcopy(ui.static_data)
    time_taken = list()

    for var_combo in new_vars:
        start = time()
        all_fits = dict()
        all_fits['yield'] = {}
        all_fits['halflife'] = {}

        csv_name = str(list(var_combo.values())).strip('[]').strip(']').replace(', ', '-') + naming_modifier

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