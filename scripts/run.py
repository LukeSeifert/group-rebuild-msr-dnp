import nlls
from radrun import Run
from simple import IrradSimple
from counts import DelayedCounts
import numpy as np
import pandas as pd
from copy import deepcopy

if __name__ == '__main__':
    import ui

    run_omc = False
    decay_daughter = True
    csv_name = 'default'

    all_fits = dict()
    all_fits['yield'] = {}
    all_fits['halflife'] = {}

    dt = 0.1
    tf = ui.default_omc_decay_time

    Count = DelayedCounts(dt, tf)
    runner = Run(ui.nuc_list,
                 run_omc=run_omc,
                 decay_track=False,
                 write_concs=False)
    dec_runner = Run(ui.nuc_list,
                 run_omc=run_omc,
                 decay_track=True,
                 write_concs=True)
    if decay_daughter:
        run_obj = dec_runner
    else:
        run_obj = runner

    data = ui.pulse_data
    a_fit, lam_fit = nlls.gen_fit(IrradSimple,
                                  data,
                                  'pulse',
                                  run_obj,
                                  Count,
                                  nlls.nlls_fit)
    all_fits['yield'][data['name']] = a_fit
    all_fits['halflife'][data['name']] = np.log(2) / lam_fit

    data = ui.static_data
    a_fit, lam_fit = nlls.gen_fit(IrradSimple,
                                  data,
                                  'saturation',
                                  run_obj,
                                  Count,
                                  nlls.nlls_fit)
    all_fits['yield'][data['name']] = a_fit
    all_fits['halflife'][data['name']] = np.log(2) / lam_fit


    base_dict = deepcopy(all_fits)
    hl_fits = base_dict['halflife']
    a_fits = base_dict['yield']
    for fit_i, fit in enumerate(a_fits.keys()):
        if fit != ui.pulse_data['name']:
            hls  = np.append(hl_fits[fit][:4], hl_fits[ui.pulse_data['name']][4:])
            ylds = np.append(a_fits[fit][:4], a_fits[ui.pulse_data['name']][4:])
            all_fits['yield'][f'{fit}-{ui.pulse_data["name"]}'] = ylds
            all_fits['halflife'][f'{fit}-{ui.pulse_data["name"]}'] = hls

    nlls.generate_csvs(all_fits, csv_name=csv_name)