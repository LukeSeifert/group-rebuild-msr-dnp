import numpy as np
from scipy.optimize import curve_fit
from scipy.optimize import least_squares
from copy import deepcopy
from counts import DelayedCounts
import os
from radrun import Run
from simple import IrradSimple
from uncertainties import ufloat, unumpy
import time
import matplotlib.pyplot as plt
import pandas as pd


class NLLS:
    """
    This class handles the non-linear least-squares functions and solvers

    """

    def __init__(self, groups: int, efficiency: float,
                 fission_term: float, times: list,
                 counts: list, a_vals_fix: list,
                 lam_vals_fix: list,
                 irradobj: IrradSimple):
        self.num_groups = groups
        self.efficiency = efficiency
        self.fission_term = fission_term
        self.times = times
        self.counts = counts
        self.a_vals_fix = a_vals_fix
        self.lam_vals_fix = lam_vals_fix
        self.fit_type = None
        self.num_unknowns_a   = len([i for i in a_vals_fix if i is None])
        self.num_unknowns_lam = len([i for i in lam_vals_fix if i is None])
        self.num_unknowns = self.num_unknowns_a + self.num_unknowns_lam
        self.t_incore = irradobj.t_incore
        self.t_excore = irradobj.t_excore
        self.t_irrad  = irradobj.net_irrad_time_s
        return
    
    def _apply_fixed_terms(self, vector_vals: list):
        a_vals = deepcopy(self.a_vals_fix)
        lam_vals = deepcopy(self.lam_vals_fix)
        vector_val_index = 0

        for ai, a in enumerate(a_vals):
            if type(a) == type(None):
                a_vals[ai] = vector_vals[vector_val_index]
                vector_val_index += 1
        
        for lami, lam in enumerate(lam_vals):
            if type(lam) == type(None):
                lam_vals[lami] = vector_vals[vector_val_index]
                vector_val_index += 1

        return a_vals, lam_vals

    
    def _group_summer(self, t, *vector_vals: list):
        group_sum = 0
        a_vals, lam_vals = self._apply_fixed_terms(vector_vals)
        for group in range(self.num_groups):
            if self.fit_type == 'pulse':
                group_val = (a_vals[group] * lam_vals[group] * np.exp(-lam_vals[group] * t))
                #group_val = (np.log(a_vals[group] * lam_vals[group]) - lam_vals[group] * t)
            elif self.fit_type == 'saturation' or self.fit_type == 'simpleflow':
                lam = lam_vals[group]
                exp_l = np.exp(-lam * t)
                exp_T = np.exp(-lam * self.t_irrad)
                exp_ex = np.exp(-lam * self.t_excore)
                sum_term = 0
                eval_time = 0
                j = 1
                while eval_time < self.t_irrad:
                    sum_term += np.exp(lam * (j * self.t_incore + (j-1) * self.t_excore - self.t_irrad))
                    j += 1
                    eval_time += self.t_excore + self.t_incore
                group_val = a_vals[group] * exp_l * (1-exp_T + (1-exp_ex) * sum_term)
                #group_val = (a_vals[group] * np.exp(-lam_vals[group] * t))
                #group_val = (np.log(a_vals[group]) - lam_vals[group] * t)
            group_sum += group_val
        #delnu = group_sum
        delnu = np.log(group_sum)
        return delnu
    
    def _group_combined_summer(self, t, *vector_vals: list):
        group_sum = 0
        vector_vals = vector_vals[0]
        a_vals = vector_vals[:6]
        lam_vals = vector_vals[6:]
        for group in range(self.num_groups):
            if self.fit_type == 'pulse':
                group_val = (a_vals[group] * lam_vals[group] * np.exp(-lam_vals[group] * t))
            elif self.fit_type == 'saturation' or self.fit_type == 'simpleflow':
                lam = lam_vals[group]
                exp_l = np.exp(-lam * t)
                exp_T = np.exp(-lam * self.t_irrad)
                exp_ex = np.exp(-lam * self.t_excore)
                #group_val = (a_vals[group] * np.exp(-lam_vals[group] * t))
                sum_term = 0
                eval_time = 0
                j = 0
                while eval_time < self.t_irrad:
                    sum_term += np.exp(lam * (j * self.t_incore + (j-1) * self.t_excore - self.t_irrad))
                    j += 1
                    eval_time += self.t_excore + self.t_incore
                group_val = a_vals * exp_l * (1-exp_T + (1-exp_ex) * sum_term)
            group_sum += group_val
        delnu = group_sum
        return delnu

    def group_fit(self, fit_type: str):
        valid_types = ['pulse', 'saturation', 'simpleflow']
        if fit_type in valid_types:
            self.fit_type = fit_type
            func = self._group_summer
        else:
            raise Exception(f'{fit_type=} not in {valid_types=}')
        
        start = time.time()
        #adjusted_counts = [i / (self.fission_term * self.efficiency) for i in self.counts]
        adjusted_counts = [unumpy.log(i / (self.fission_term * self.efficiency)) for i in self.counts]
        adjusted_count_vals = [unumpy.nominal_values(x) for x in adjusted_counts]
        adjusted_count_uncerts = [unumpy.std_devs(x) for x in adjusted_counts]
        
        #params, covariance, info, _, _ = curve_fit(func, self.times, adjusted_counts,
        #                            p0=[1]*self.num_unknowns, maxfev=100000,
        #                            bounds=(0, 1e3), full_output=True,
        #                            xtol=2.23e-16, gtol=2.23e-16,
        #                            verbose=0, ftol=2.23e-16)
        params, covariance, info, _, _ = curve_fit(func, self.times, adjusted_count_vals,
                                    p0=[1]*self.num_unknowns,
                                    method='trf',
                                    bounds=(0, 1e3), full_output=True,
                                    maxfev=1e6,
                                    #sigma=adjusted_count_uncerts,
                                    gtol=None,
                                    xtol=None,
                                    verbose=0,
                                    ftol=2.23e-16)
        end = time.time()
        print(f'Took {round(end-start, 3)}s for NLLS fit')

        chi_squared = np.sum(info['fvec']**2)
        print(f'{chi_squared=}')

        a_fits = self.a_vals_fix
        a_vals = params[:self.num_unknowns_a].tolist()
        a_counter = 0
        for ai, a in enumerate(a_fits):
            if a == None:
                a_fits[ai] = a_vals[a_counter]
                a_counter += 1

        lam_fits = self.lam_vals_fix
        lam_vals = params[self.num_unknowns_a:].tolist()
        lam_counter = 0
        for lami, lam in enumerate(lam_fits):
            if lam == None:
                lam_fits[lami] = lam_vals[lam_counter]
                lam_counter += 1

        zipped = list(zip(a_fits, lam_fits))
        sorted_zipped = sorted(zipped, key=lambda x: x[1])
        a_fits_sorted, lam_fits_sorted = zip(*sorted_zipped)
        a_fits = list(a_fits_sorted)
        lam_fits = list(lam_fits_sorted)

        print(f'{np.linalg.cond(covariance)=}')
        self.fit_type = None
        return a_fits, lam_fits
    
    def _plot(self, name: str, times: list, counts: list,
              a_fits: list, lam_fits: list, fit_type: str):
        params = a_fits + lam_fits
        group_counts = list()
        self.fit_type = fit_type
        for t in times:
            #delnu = (self.fission_term * self.efficiency * 
            #         self._group_summer(t, params))
            delnu = (self.fission_term * self.efficiency * 
                     np.exp(self._group_summer(t, params)))
            group_counts.append(delnu)
        use_counts = np.asarray([unumpy.nominal_values(x) for x in counts])
        d_counts = np.asarray([unumpy.std_devs(x) for x in counts])
        plt.plot(times, use_counts, label='Count data')
        plt.fill_between(times, use_counts+d_counts, use_counts-d_counts, alpha=0.6)
        plt.plot(times, group_counts, label='Group data')
        plt.xlabel('Time [s]')
        plt.yscale('log')
        plt.ylabel('Delayed Neutron Count Rate')
        plt.legend()
        plt.tight_layout()
        save_path = f'./images/{name}_groupcompare.png'
        try:
            plt.savefig(save_path)
        except FileNotFoundError:
            os.mkdir('./images')
            plt.savefig(save_path)
        plt.close()

        pcnt_diff = [(use_counts[i] - group_counts[i])/(use_counts[i])* 100 for i in range(len(times))]
        plt.plot(times, pcnt_diff)
        plt.xlabel('Time [s]')
        plt.ylabel('Percent Difference')
        plt.tight_layout()
        save_path = f'./images/{name}_pcntdiff.png'
        plt.savefig(save_path)
        plt.close()


        self.fit_type = None
        return



def _print_helper(name, a_fits, tot_yield, lam_fits, half_lives):
    print(f'{name=}')
    print(f'    yields = {np.round(a_fits, 5).tolist()}')
    print(f'    hls = {np.round(half_lives, 5).tolist()}')
    print(f'    lams = {np.round(lam_fits, 5).tolist()}')
    print(f'    delnu (from summed yields) = {float(np.round(tot_yield, 5))}')
    print()
    return


def keepin_test(Count: DelayedCounts):
    name = 'Keepin Fit (Pulse)'
    yields = [0.00063, 0.00351, 0.00310, 0.00672, 0.00211, 0.00043]
    hls = [54.51, 21.84, 6.00, 2.23, 0.496, 0.179]
    lams = [np.log(2)/hl for hl in hls]
    fissions = 1e16
    times, counts = Count.from_groups(yields, lams, fissions)
    
    num_groups = 6
    a_vals_fix = [None] * 6
    lam_vals_fix = [None] * 6
    group = NLLS(groups=num_groups, efficiency=1, fission_term=fissions,
                times=times, counts=counts, a_vals_fix=a_vals_fix,
                lam_vals_fix=lam_vals_fix)
    a_fits, lam_fits = group.group_fit('pulse')
    half_lives = [np.log(2)/lam for lam in lam_fits]
    tot_yield = sum(a_fits)
    _print_helper(name, a_fits, tot_yield, lam_fits, half_lives)
    return a_fits, lam_fits

def from_counts(name: str, fission_term: float, Count: DelayedCounts,
                a_vals_fix: list, lam_vals_fix: list,
                irrad_type: str,
                output_path: str,
                irradobj: IrradSimple,
                cutoff_scale: float=1):
    num_groups = len(a_vals_fix)
    csv_path = f'{output_path}/concs.csv'
    times, counts = Count.from_concs(csv_path, cutoff_scale=cutoff_scale)

    num_groups = 6
    group = NLLS(groups=num_groups, efficiency=1, fission_term=fission_term,
                 times=times, counts=counts, a_vals_fix=a_vals_fix,
                 lam_vals_fix=lam_vals_fix, irradobj=irradobj)
    a_fits, lam_fits = group.group_fit(irrad_type)
    group._plot(name, times, counts, a_fits, lam_fits, irrad_type)
    half_lives = [np.log(2)/lam for lam in lam_fits]
    tot_yield = sum(a_fits)
    _print_helper(name, a_fits, tot_yield, lam_fits, half_lives)

    return a_fits, lam_fits

def group_combined_fit(groups):
    pulse_func = groups[0]._group_combined_summer
    pulse_times = groups[0].times
    pulse_counts = groups[0].counts
    groups[0].fit_type = 'pulse'
    pulse_counts = [i.n / (groups[0].fission_term * groups[0].efficiency) for i in groups[0].counts]

    sat_func = groups[1]._group_combined_summer
    sat_times = groups[1].times
    sat_counts = groups[1].counts
    sat_counts = [i.n / (groups[1].fission_term * groups[1].efficiency) for i in groups[1].counts]
    groups[1].fit_type = 'saturation'

    p0=[1]*groups[0].num_unknowns
    #adjusted_counts = [unumpy.log(i / (self.fission_term * self.efficiency)) for i in self.counts]
    #adjusted_count_vals = [unumpy.nominal_values(x) for x in adjusted_counts]

    def residual_func(parameters, pulse_times, pulse_counts, sat_times, sat_counts):
        parameters = np.array(parameters, dtype=float)
        pulse_residual = (pulse_counts - pulse_func(pulse_times, parameters)) / pulse_counts
        sat_residual = (sat_counts - sat_func(sat_times, parameters)) / sat_counts
        net_residual = pulse_residual + sat_residual
        print(net_residual)
        return net_residual

    start = time.time()
    result = least_squares(residual_func, p0, bounds=(0, 1000), method='trf',
                  ftol=None, xtol=None, gtol=1e-8,
                  verbose=2,
                  args=(pulse_times, pulse_counts, sat_times, sat_counts))
    print(f'Took {round(time.time()-start, 3)}s for NLLS fit')
    params = result.x
    print(result)
    print(f'{params = }')
    #print(residual_func(result.x, pulse_times, pulse_counts, sat_times, sat_counts))
    


    a_fits = groups[0].a_vals_fix
    a_vals = params[:groups[0].num_unknowns_a].tolist()
    a_counter = 0
    for ai, a in enumerate(a_fits):
        if a == None:
            a_fits[ai] = a_vals[a_counter]
            a_counter += 1

    lam_fits = groups[0].lam_vals_fix
    lam_vals = params[groups[0].num_unknowns_a:].tolist()
    lam_counter = 0
    for lami, lam in enumerate(lam_fits):
        if lam == None:
            lam_fits[lami] = lam_vals[lam_counter]
            lam_counter += 1

    zipped = list(zip(a_fits, lam_fits))
    sorted_zipped = sorted(zipped, key=lambda x: x[1])
    a_fits_sorted, lam_fits_sorted = zip(*sorted_zipped)
    a_fits = list(a_fits_sorted)
    lam_fits = list(lam_fits_sorted)
    print(sum(a_fits))

    return a_fits, lam_fits
    


def from_combined_counts(names: list, fission_terms: list, Counts: list,
                a_vals_fix: list, lam_vals_fix: list,
                irrad_types: list,
                output_paths: list,
                irrad_objs: list,
                cutoff_scale: float=1):
    num_groups = len(a_vals_fix)
    groups = list()
    for i in range(len(names)):
        csv_path = f'{output_paths[i]}/concs.csv'
        times, counts = Counts[i].from_concs(csv_path, cutoff_scale=cutoff_scale)

        group = NLLS(groups=num_groups, efficiency=1, fission_term=fission_terms[i],
                    times=times, counts=counts, a_vals_fix=a_vals_fix,
                    lam_vals_fix=lam_vals_fix, irradobj=irrad_objs[i])
        groups.append(group)
    a_fits, lam_fits = group_combined_fit(groups)

    return a_fits, lam_fits

def nlls_fit(IrradObj: IrradSimple, irrad_type: str, runner: Run,
             Count: DelayedCounts, num_groups=6):
    name = IrradObj.name
    output_path = IrradObj.output_path
    avgF, netF = runner.simple_compare(IrradObj)
    runner._reset_metadict()
    if irrad_type == 'pulse':
        fission_term = netF
    elif irrad_type == 'simpleflow' or irrad_type == 'saturation':
        fission_term = netF / IrradObj.net_irrad_time_s #avgF
        #fission_term = avgF

    #yields = [0.0004, 0.00171, 0.00245, 0.00075]
    #lams = np.log(2) / [12.71955, 4.64925, 1.86579, 0.34046]
    #a_vals_fix = [None] * 2 + yields
    #lam_vals_fix = [None] * 2 + list(lams)

    a_vals_fix = [None] * num_groups
    lam_vals_fix = [None] * num_groups
    cutoff_scale = 1
    a_fits, lam_fits = from_counts(name, fission_term, Count,
                                   a_vals_fix, lam_vals_fix,
                                   irrad_type,
                                   output_path,
                                   IrradObj,
                                   cutoff_scale)
    return a_fits, lam_fits


def nlls_combined_fit(IrradObj: list, irrad_type: list, runner: list,
                      Count: list, num_groups=6):
    names = list()
    output_paths = list()
    fission_terms = list()
    for i, irradobj in enumerate(IrradObj):
        name = irradobj.name
        output_path = irradobj.output_path
        avgF, netF = runner[i].simple_compare(irradobj)
        runner[i]._reset_metadict()
        if irrad_type[i] == 'pulse':
            fission_term = netF
        elif irrad_type[i] == 'simpleflow' or irrad_type[i] == 'saturation':
            fission_term = netF / irradobj.net_irrad_time_s #avgF
        else:
            raise ValueError(f'{irrad_type[i]} invalid')
        names.append(name)
        output_paths.append(output_path)
        fission_terms.append(fission_term)

    a_vals_fix = [None] * num_groups
    lam_vals_fix = [None] * num_groups
    cutoff_scale = 1
    a_fits, lam_fits = from_combined_counts(names, fission_terms, Count,
                                   a_vals_fix, lam_vals_fix,
                                   irrad_type,
                                   output_paths,
                                   IrradObj,
                                   cutoff_scale)
    return a_fits, lam_fits

def gen_fit(irrad_obj : IrradSimple,
            irrad_type : str,
            runner : Run,
            Count : DelayedCounts,
            fit_func):

    a_fits, lam_fits = fit_func(irrad_obj, irrad_type, runner, Count)
    return a_fits, lam_fits

def generate_csvs(all_fits: dict,
                  csv_path: str = './postprocess',
                  csv_name: str = 'default'):
    yld_data = {'Yield': [],
                'Group': [],
                'Data Source': []}
    hls_data = {'Half-life [s]': [],
                'Group': [],
                'Data Source': []}
    for fit_type in all_fits.keys():
        for source in all_fits[fit_type].keys():
            use_data = all_fits[fit_type][source]
            for group in range(len(use_data)):
                if fit_type == 'yield':
                    yld_data['Yield'].append(use_data[group])
                    yld_data['Group'].append(group+1)
                    yld_data['Data Source'].append(source)
                elif fit_type == 'halflife':
                    hls_data['Half-life [s]'].append(use_data[group])
                    hls_data['Group'].append(group+1)
                    hls_data['Data Source'].append(source)
                else:
                    raise KeyError(f'Key {fit_type} not available as option')
    df = pd.DataFrame(yld_data)
    df.to_csv(f'{csv_path}/yields/{csv_name}.csv')

    df = pd.DataFrame(hls_data)
    df.to_csv(f'{csv_path}/halflives/{csv_name}.csv')
    return

if __name__ == "__main__":
    import ui
    run_omc = True


    dt = 0.1
    tf = ui.default_final_time

    Count = DelayedCounts(dt, tf)
    runner = Run(ui.nuc_list,
                 run_omc=run_omc,
                 decay_track=False,
                 write_concs=False)
    dec_runner = Run(ui.nuc_list,
                 run_omc=run_omc,
                 decay_track=True,
                 write_concs=True)

#    a_fits, lam_fits = keepin_test(Count)

    irradobj = IrradSimple(data_dict=ui.pulse_data)
    a_fits, lam_fits = nlls_fit(irradobj, 'pulse', dec_runner, Count)

#    irradobj = IrradSimple(data_dict=ui.static_data)
#    a_fits, lam_fits = nlls_fit(irradobj, 'saturation', runner, Count)
#    a_fits, lam_fits = nlls_fit(irradobj, 'saturation', dec_runner, Count)

#    irradobj = IrradSimple(data_dict=ui.flow_repr_data)
#    a_fits, lam_fits = nlls_fit(irradobj, 'simpleflow', runner, Count)

#    irradobj = IrradSimple(data_dict=ui.flow_data)
#    a_fits, lam_fits = nlls_fit(irradobj, 'simpleflow', runner, Count)
#
#    irradobj = IrradSimple(data_dict=ui.mostly_excore_data)
#    a_fits, lam_fits = nlls_fit(irradobj, 'simpleflow', runner, Count)
#
#
#    irradobj = IrradSimple(data_dict=ui.exflow_repr_data)
#    a_fits, lam_fits = nlls_fit(irradobj, 'simpleflow', runner, Count)
#
#    irradobj = IrradSimple(data_dict=ui.static_repr_data)
#    a_fits, lam_fits = nlls_fit(irradobj, 'saturation', runner, Count)
