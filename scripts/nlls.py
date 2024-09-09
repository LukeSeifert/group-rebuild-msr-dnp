import numpy as np
from scipy.optimize import curve_fit

def gaussian(x, a, b, c):
    return a * np.exp(-((x - b) ** 2) / (2 * c ** 2))

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

        return
    
    def _apply_fixed_terms(self, a_vals, lam_vals):
        for ai, a in enumerate(self.a_vals_fix):
            if type(a) != type(None):
                a_vals[ai] = a
        
        for lami, lam in enumerate(self.lam_vals_fix):
            if type(lam) != type(None):
                lam_vals[lami] = lam
        
        return a_vals, lam_vals
    
    def _pulse(self, t, a_vals, lam_va        a_vals, lam_vals = self._apply_fixed_terms(a_vals, lam_vals)ls):
        group_sum = 0
        a_vals, lam_vals = self._apply_fixed_terms(a_vals, lam_vals)
        for group in range(self.num_groups):
            group_sum += (a_vals[group] * lam_vals[group] * 
                          np.exp(-lam_vals[group] * t))
        delnu = self.efficiency * self.fission_term * group_sum
        return delnu
    
    def _saturation(self, t, a_vals, lam_vals):
        group_sum = 0
        a_vals, lam_vals = self._apply_fixed_terms(a_vals, lam_vals)
        for group in range(self.num_groups):
            group_sum += (a_vals[group] * np.exp(-lam_vals[group] * t))
        delnu = self.efficiency * self.fission_term * group_sum
        return delnu
    
    def curve_fit(self, fit_type: str):
        if fit_type == 'pulse':
            func = self._pulse
        elif fit_type == '_saturation':
            func = self._saturation
        else:
            raise Exception(f'{fit_type=} not implemented')
        
        params, covariance = curve_fit(func, self.times, self.counts)
        a_fits = params[:self.num_groups]
        lam_fits = params[self.num_groups:]
        return a_fits, lam_fits


if __name__ == "__main__":
    import 