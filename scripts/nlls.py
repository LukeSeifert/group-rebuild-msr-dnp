import numpy as np
from scipy.optimize import curve_fit
from copy import deepcopy
from counts import DelayedCounts
import os


class NNLS:
    """
    This class handles the non-linear least-squares functions and solvers

    """

    def __init__(self, groups: int, efficiency: float,
                 fission_term: float, times: list,
                 counts: list, a_vals_fix: list,
                 lam_vals_fix: list, linearize: bool):
        self.num_groups = groups
        self.efficiency = efficiency
        self.fission_term = fission_term
        self.times = times
        self.counts = counts
        self.a_vals_fix = a_vals_fix
        self.lam_vals_fix = lam_vals_fix
        self.linearize = linearize
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
            elif self.fit_type == 'saturation' or self.fit_type == 'simpleflow':
                group_val = (a_vals[group] * np.exp(-lam_vals[group] * t))
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
        
        adjusted_counts = [i / (self.fission_term * self.efficiency) for i in self.counts]
        
        params, covariance, info, _, _ = curve_fit(func, self.times, adjusted_counts,
                                    p0=[1]*self.num_unknowns, maxfev=100000,
                                    bounds=(0, 1e3), full_output=True,
                                    xtol=None, gtol=None,
                                    verbose=0, sigma=np.asarray(adjusted_counts))

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
            delnu = (self.fission_term * self.efficiency * 
                     self._group_summer(t, params))
            group_counts.append(delnu)
        plt.plot(times, counts, label='Count data')
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

        pcnt_diff = [(counts[i] - group_counts[i])/(counts[i])* 100 for i in range(len(times))]
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
    print(f'    delnu = {float(np.round(tot_yield, 5))}')
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
                lam_vals_fix=lam_vals_fix, linearize=False)
    a_fits, lam_fits = group.group_fit('pulse')
    half_lives = [np.log(2)/lam for lam in lam_fits]
    tot_yield = sum(a_fits)
    _print_helper(name, a_fits, tot_yield, lam_fits, half_lives)
    return a_fits, lam_fits

def from_counts(name: str, fission_term: float, Count: DelayedCounts,
                a_vals_fix: list, lam_vals_fix: list,
                irrad_type: str,
                linearize: bool,
                cutoff_scale: float=1):
    num_groups = len(a_vals_fix)
    csv_path = f'./results/{name}/concs.csv'
    times, counts = Count.from_concs(csv_path, cutoff_scale=cutoff_scale)

    num_groups = 6
    group = NNLS(groups=num_groups, efficiency=1, fission_term=fission_term,
                 times=times, counts=counts, a_vals_fix=a_vals_fix,
                 lam_vals_fix=lam_vals_fix, linearize=linearize)
    a_fits, lam_fits = group.group_fit(irrad_type)
    group._plot(name, times, counts, a_fits, lam_fits, irrad_type)
    half_lives = [np.log(2)/lam for lam in lam_fits]
    tot_yield = sum(a_fits)
    _print_helper(name, a_fits, tot_yield, lam_fits, half_lives)


    return a_fits, lam_fits



if __name__ == "__main__":
    import matplotlib.pyplot as plt
    import ui


    dt = 1e-1
    tf = 1000
    Count = DelayedCounts(dt, tf)
    linearize = True

    #a_fits, lam_fits = keepin_test(Count)

    
    name = 'Pulse'
    irrad_type = 'pulse'
    net_fiss = 8.502E+14
    num_groups = 6
    fission_term = net_fiss
    a_vals_fix = [None] * num_groups
    lam_vals_fix = [None] * num_groups
    cutoff_scale = 1
    pulse_a_fits, pulse_lam_fits = from_counts(name, fission_term, Count,
                                               a_vals_fix, lam_vals_fix,
                                               irrad_type, linearize,
                                               cutoff_scale)
    
    Count = DelayedCounts(dt, tf, t0=0)
    name = 'Static'
    irrad_type = 'saturation'
    avg_fiss_rate = 4.290E+13
    num_groups = 6
    fission_term = avg_fiss_rate
    #a_vals_fix = [None] * 2 + pulse_a_fits[2:]
    #lam_vals_fix = [None] * 2 + pulse_lam_fits[2:]
    a_vals_fix = [None] * 6
    lam_vals_fix = [None] * 6
    cutoff_scale = 1
    a_fits, lam_fits = from_counts(name, fission_term, Count, a_vals_fix, 
                                   lam_vals_fix, irrad_type, linearize,
                                   cutoff_scale)
    

    name = 'Flowing'
    irrad_type = 'simpleflow'
    avg_fiss_rate = 4.222E+13
    num_groups = 6
    fission_term = avg_fiss_rate
    a_vals_fix = [None] * 2 + pulse_a_fits[2:]
    lam_vals_fix = [None] * 2 + pulse_a_fits[2:]
    cutoff_scale = 1
    a_fits, lam_fits = from_counts(name, fission_term, Count, a_vals_fix, 
                                   lam_vals_fix, irrad_type, linearize,
                                   cutoff_scale)



