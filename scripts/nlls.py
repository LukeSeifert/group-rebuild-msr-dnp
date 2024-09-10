import numpy as np
from scipy.optimize import curve_fit


class NNLS:
    """
    This class handles the non-linear least-squares functions and solvers

    """

    def __init__(self, groups: int, efficiency: float,
                 fission_term: float, times: list,
                 counts: list, a_vals_fix: list,
                 lam_vals_fix: list):
        # TODO - allow some groups to be hardcoded
        self.num_groups = groups
        self.efficiency = efficiency
        self.fission_term = fission_term
        self.times = times
        self.counts = counts
        self.a_vals_fix = a_vals_fix
        self.lam_vals_fix = lam_vals_fix
        self.fit_type = None

        return
    
    def _apply_fixed_terms(self, vector_vals: list):
        a_vals = np.zeros(self.num_groups)
        lam_vals = np.zeros(self.num_groups)

        for ai in range(self.num_groups):
            a = self.a_vals_fix[ai]
            if type(a) != type(None):
                a_vals[ai] = a
            else:
                a_vals[ai] = vector_vals[ai]
        
        for lami in range(self.num_groups, self.num_groups*2):
            lam = self.lam_vals_fix[lami-self.num_groups]
            if type(lam) != type(None):
                lam_vals[lami] = lam
            else:
                lam_vals[lami-self.num_groups] = vector_vals[lami]
        
        return a_vals, lam_vals

    
    def _group_summer(self, t, *vector_vals: list,
                      kappa:int = 1, simpleflow: bool=False):
        group_sum = 0
        a_vals, lam_vals = self._apply_fixed_terms(vector_vals)
        for group in range(self.num_groups):
            if self.fit_type == 'pulse':
                group_val = (a_vals[group] * lam_vals[group] * np.exp(-lam_vals[group] * t))
            elif self.fit_type == 'saturation' or self.fit_type == 'simpleflow':
                group_val = (a_vals[group] * np.exp(-lam_vals[group] * t))
            eta_sum = 0
            eta_sum = self._get_eta_sum(kappa, simpleflow, a_vals[group],
                                          lam_vals[group])
            group_sum += group_val * eta_sum
        delnu = self.efficiency * self.fission_term * group_sum
        return delnu

    def group_fit(self, fit_type: str):
        valid_types = ['pulse', 'saturation', 'simpleflow']
        if fit_type in valid_types:
            self.fit_type = fit_type
            func = self._group_summer
            self.fit_type = None
        else:
            raise Exception(f'{fit_type=} not in {valid_types=}')
        
        params, covariance = curve_fit(func, self.times, self.counts,
                                       p0=[1]*self.num_groups*2, maxfev=100000,
                                       bounds=(0, 1e3))
        a_fits = params[:self.num_groups]
        lam_fits = params[self.num_groups:]

        zipped = list(zip(a_fits, lam_fits))
        sorted_zipped = sorted(zipped, key=lambda x: x[1])
        a_fits_sorted, lam_fits_sorted = zip(*sorted_zipped)
        a_fits = list(a_fits_sorted)
        lam_fits = list(lam_fits_sorted)

        print(f'{np.linalg.cond(covariance)=}')
        return a_fits, lam_fits

if __name__ == "__main__":
    import matplotlib.pyplot as plt
    from counts import DelayedCounts
    import ui
    dt = 1e-1
    tf = 1000
    cutoff = 1


    time_list = list()
    count_list = list()
    name_list = list()
    fission_list = list()
    Count = DelayedCounts(dt, tf)


    name = 'Keepin Fit (Pulse)'
    yields = [0.00063, 0.00351, 0.00310, 0.00672, 0.00211, 0.00043]
    hls = [54.51, 21.84, 6.00, 2.23, 0.496, 0.179]
    lams = [np.log(2)/hl for hl in hls]
    fissions = 1E16
    times, counts = Count.from_groups(yields, lams, fissions)
    
    num_groups = 6
    group = NNLS(groups=num_groups, efficiency=1, fission_term=fissions,
                 times=times, counts=counts, a_vals_fix=[None]*num_groups,
                 lam_vals_fix=[None]*num_groups)
    a_fits, lam_fits = group.group_fit('pulse')
    half_lives = [np.log(2)/lam for lam in lam_fits]
    print(f'{np.round(a_fits, 5)  =}')
    print(f'{np.round(lam_fits, 5)=}')
    print(f'{np.round(half_lives, 5)=}')


