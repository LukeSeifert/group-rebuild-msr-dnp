import numpy as np
from scipy.optimize import curve_fit
from copy import deepcopy
from counts import DelayedCounts
import os
from radrun import Run
from simple import IrradSimple
from uncertainties import ufloat, unumpy
import time


class NNLS:
    """
    This class handles the non-linear least-squares functions and solvers

    """

    def __init__(self, groups: int, efficiency: float,
                 fission_term: float, times: list,
                 counts: list, a_vals_fix: list,
                 lam_vals_fix: list):
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
                group_val = (a_vals[group] * np.exp(-lam_vals[group] * t))
                #group_val = (np.log(a_vals[group]) - lam_vals[group] * t)
            group_sum += group_val
        #delnu = group_sum
        delnu = np.log(group_sum)
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
                                    p0=[1]*self.num_unknowns, maxfev=100000,
                                    bounds=(0, 1e3), full_output=True,
                                    #sigma=adjusted_count_uncerts,
                                    xtol=2.23e-16, gtol=2.23e-16,
                                    verbose=0, ftol=2.23e-16)
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
    group = NNLS(groups=num_groups, efficiency=1, fission_term=fissions,
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
                cutoff_scale: float=1):
    num_groups = len(a_vals_fix)
    csv_path = f'./results/{name}/concs.csv'
    times, counts = Count.from_concs(csv_path, cutoff_scale=cutoff_scale)

    num_groups = 6
    group = NNLS(groups=num_groups, efficiency=1, fission_term=fission_term,
                 times=times, counts=counts, a_vals_fix=a_vals_fix,
                 lam_vals_fix=lam_vals_fix)
    a_fits, lam_fits = group.group_fit(irrad_type)
    group._plot(name, times, counts, a_fits, lam_fits, irrad_type)
    half_lives = [np.log(2)/lam for lam in lam_fits]
    tot_yield = sum(a_fits)
    _print_helper(name, a_fits, tot_yield, lam_fits, half_lives)


    return a_fits, lam_fits

def nlls_fit(IrradObj: IrradSimple, irrad_type: str, runner: Run,
             Count: DelayedCounts):
    name = IrradObj.name
    avgF, netF = runner.simple_compare(IrradObj)
    runner._reset_metadict()
    num_groups = 6
    if irrad_type == 'pulse':
        fission_term = netF
    elif irrad_type == 'simpleflow' or irrad_type == 'saturation':
        fission_term = avgF

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
                                   cutoff_scale)
    return a_fits, lam_fits


if __name__ == "__main__":
    import matplotlib.pyplot as plt
    import ui


    dt = 1e-1
    tf = 500
    Count = DelayedCounts(dt, tf)
    runner = Run(ui.nuc_list,
                 run_omc=False,
                 decay_track=True,
                 write_concs=True)

#    a_fits, lam_fits = keepin_test(Count)

#    irradobj = IrradSimple(data_dict=ui.pulse_data)
#    a_fits, lam_fits = nlls_fit(irradobj, 'pulse', runner, Count)

    ui.static_data['name'] = ui.static_data['name'] #+ 'Decay'
    irradobj = IrradSimple(data_dict=ui.static_data)
    a_fits, lam_fits = nlls_fit(irradobj, 'saturation', runner, Count)

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