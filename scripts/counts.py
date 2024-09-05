import numpy as np
from simple import IrradSimple
import pandas as pd
import matplotlib.pyplot as plt
import os
import time
import scipy as sp


class DelayedCounts:
    """
    This class generates delayed neutron counts for given times either from
        concentrations and IAEA data, or from provided group data    
    """
    def __init__(self, dt: float, tf: float):
        self.dt = dt
        self.tf = tf
        self.times = np.arange(0, tf+dt, dt)

        simple_irrad_obj = IrradSimple(None)
        self.iaea_data = simple_irrad_obj._get_data()
        return
    
    def _order_iaea_data(self):
        self.pns = dict()
        self.lams = dict()
        for i, nuc in enumerate(self.iaea_data['nucid']):
            beta_prob = self.iaea_data['  beta- %'][i]
            Pn1 = self.iaea_data[' pn1 % '][i]
            Pn2 = self.iaea_data[' pn2 % '][i]
            Pn3 = self.iaea_data[' pn3 % '][i]
            prob_neutron = beta_prob * (1*Pn1 + 2*Pn2 + 3*Pn3)
            lam = self.iaea_data[' T1/2 [s] '][i]
            self.pns[nuc] = prob_neutron
            self.lams[nuc] = np.log(2) / lam 
        return
    
    def from_concs(self, csv_path, cutoff_scale = 1):
        start = time.time()
        self._order_iaea_data()
        df = pd.read_csv(csv_path)
        omc_nucs = df.iloc[:, 0]
        concs = df.iloc[:, 1:]
        use_nucs = list()
        if concs.shape[1] > 1:
            raise Exception('Concentrations as a function of time not implemented')
        counts = list()
        mult_term = dict()
        for nuc in self.iaea_data['nucid']:
            if nuc in list(omc_nucs):
                nuc_index = df[df.iloc[:, 0] == nuc].index[0]
                conc = df.iloc[nuc_index, 1]
                mult_val = self.pns[nuc] * self.lams[nuc] * conc
                if mult_val < cutoff_scale:
                    continue
                mult_term[nuc] = mult_val
                use_nucs.append(nuc)

        for t in self.times:
            run_count = 0
            for nuc in use_nucs:
                run_count += mult_term[nuc] * np.exp(-self.lams[nuc] * t)
            counts.append(run_count)

        print(f'Number of nuclides: {len(use_nucs)}')
        end = time.time()
        print(f'Took {round(end-start, 3)}s for {len(self.times)} steps')
        return self.times, counts
    
    def from_groups(self, yields: list, lams: list, fissions: float,
                    irradiation_time: float=None):
        counts = list()
        if type(irradiation_time) is type(None):
            leading_term = fissions
            mult_term = [yields[i] * lams[i] for i in range(len(yields))]
        else:
            leading_term = fissions / irradiation_time
            mult_term = yields
        for t in self.times:
            run_count = 0
            for group in range(len(yields)):
                run_count += mult_term[group] * np.exp(-lams[group] * t)
            counts.append(leading_term * run_count)
        print(f'Group delnu: {round(sum(lams), 3)}')
        return self.times, counts
    
    def get_delnu(self, fissions: float, counts: list, times: list):
        """
        Generate delayed neutrons per fission. Only for Pulse irradiation.

        Parameters
        ----------
        fissions : float
            Number of fissions in pulse
        counts : list
            Counts over time
        times : list
            Times

        Returns
        -------
        delnu : float
            Delayed neutrons per fission
        """
        delnu = sp.integrate.simpson(counts, x=times) / fissions
        return delnu
        

    
    def count_compare(self, count_names: list, times: list, counts: list,
                      fissions: list = None):
        for i in range(len(count_names)):
            plt.plot(times[i], counts[i], label=count_names[i])
        plt.legend()
        plt.xlabel('Time [s]')
        plt.ylabel('Delayed Neutron Counts')
        plt.yscale('log')
        plt.tight_layout()
        try:
            plt.savefig(f'./images/counts.png')
        except FileNotFoundError:
            os.mkdir('./images')
            plt.savefig(f'./images/counts.png')
        plt.close()

        if type(fissions) != type(None):
            for i in range(len(count_names)):
                use_counts = [counts[i][j]/fissions[i] for j in range(len(counts[i]))]
                plt.plot(times[i], use_counts, label=count_names[i])
            plt.legend()
            plt.xlabel('Time [s]')
            plt.ylabel('Delayed Neutron Counts per Fission')
            plt.yscale('log')
            plt.tight_layout()
            plt.savefig(f'./images/counts_per_fiss.png')
            plt.close()
        return





if __name__ == "__main__":
    import ui
    dt = 1e-1
    tf = 1000
    cutoff = 1


    time_list = list()
    count_list = list()
    name_list = list()
    fission_list = list()
    Count = DelayedCounts(dt, tf)

    name = 'static'
    csv_path = f'./results/{name}/concs.csv'
    fissions = 1.287E12
    times, delnu_counts = Count.from_concs(csv_path, cutoff_scale=cutoff)
    delnu = Count.get_delnu(fissions, delnu_counts, times)
    print(f'{name=} {delnu=:.3E}')
    time_list.append(times)
    count_list.append(delnu_counts)
    name_list.append(name)
    fission_list.append(fissions)

    name = 'flowing'
    csv_path = f'./results/{name}/concs.csv'
    fissions = 8.318E11
    times, delnu_counts = Count.from_concs(csv_path, cutoff_scale=cutoff)
    delnu = Count.get_delnu(fissions, delnu_counts, times)
    print(f'{name=} {delnu=:.3E}')
    time_list.append(times)
    count_list.append(delnu_counts)
    name_list.append(name)
    fission_list.append(fissions)

    name = 'pulse'
    csv_path = f'./results/{name}/concs.csv'
    fissions = 8.610E10
    times, delnu_counts = Count.from_concs(csv_path, cutoff_scale=cutoff)
    delnu = Count.get_delnu(fissions, delnu_counts, times)
    print(f'{name=} {delnu=:.3E}')
    time_list.append(times)
    count_list.append(delnu_counts)
    name_list.append(name)
    fission_list.append(fissions)

    name = 'Keepin Fit (Pulse)'
    yields = [0.00063, 0.00351, 0.00310, 0.00672, 0.00211, 0.00043]
    hls = [54.51, 21.84, 6.00, 2.23, 0.496, 0.179]
    lams = [np.log(2)/hl for hl in hls]
    fissions = 8.610E10 # pulse
    times, delnu_counts = Count.from_groups(yields, lams, fissions)
    delnu = Count.get_delnu(fissions, delnu_counts, times)
    print(f'{name=} {delnu=:.3E}')
    time_list.append(times)
    count_list.append(delnu_counts)
    name_list.append(name)
    fission_list.append(fissions)

    name = 'Keepin Fit (Saturation)'
    yields = [0.00063, 0.00351, 0.00310, 0.00672, 0.00211, 0.00043]
    hls = [54.51, 21.84, 6.00, 2.23, 0.496, 0.179]
    lams = [np.log(2)/hl for hl in hls]
    fissions = 1.287E12 # static
    rad_time = 5 * 60
    times, delnu_counts = Count.from_groups(yields, lams, fissions,
                                            irradiation_time=rad_time)
    delnu = Count.get_delnu(fissions, delnu_counts, times)
    print(f'{name=} {delnu=:.3E}')
    time_list.append(times)
    count_list.append(delnu_counts)
    name_list.append(name)
    fission_list.append(fissions)
    


    Count.count_compare(count_names=name_list,
                        times=time_list,
                        counts=count_list,
                        fissions=fission_list)

