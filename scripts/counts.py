import numpy as np
from simple import IrradSimple
import pandas as pd
import matplotlib.pyplot as plt
import os
import time
import scipy as sp
from uncertainties import ufloat, unumpy


class DelayedCounts:
    """
    This class generates delayed neutron counts for given times either from
        concentrations and IAEA data, or from provided group data    
    """
    def __init__(self, dt: float, tf: float, t0: float=0):
        self.dt = dt
        self.tf = tf
        self.times = np.arange(t0, tf+dt, dt)

        simple_irrad_obj = IrradSimple(None)
        self.iaea_data = simple_irrad_obj._get_data()
        return
    
    def _order_iaea_data(self):
        self.pns = dict()
        self.lams = dict()
        for i, nuc in enumerate(self.iaea_data['nucid']):
            beta_prob = self.iaea_data['  beta- %'][i] / 100
            d_beta_prob = self.iaea_data[' D beta-'][i] / 100
            Pn1 = self.iaea_data[' pn1 % '][i] / 100
            d_Pn1 = self.iaea_data['D pn1'][i] / 100
            Pn2 = self.iaea_data[' pn2 % '][i] / 100
            d_Pn2 = self.iaea_data['D pn2 '][i] / 100
            Pn3 = self.iaea_data[' pn3 % '][i] / 100
            d_Pn3 = self.iaea_data['D pn3'][i] / 100
            prob_neutron = beta_prob * (1*Pn1 + 2*Pn2 + 3*Pn3)
            d_prob_neutron = np.sqrt(((Pn1+2*Pn2+3*Pn3) * d_beta_prob)**2 + 
                                     ((beta_prob) * d_Pn1)**2 + 
                                     ((2 * beta_prob) * d_Pn2)**2 + 
                                     ((3*beta_prob) * d_Pn3)**2)
            hl = self.iaea_data[' T1/2 [s] '][i]
            d_hl = self.iaea_data[' D T1/2 [s]'][i]
            self.pns[nuc] = ufloat(prob_neutron, d_prob_neutron)
            hl = ufloat(hl, d_hl)
            lam = unumpy.log(2) / hl
            self.lams[nuc] = lam
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
                run_count += mult_term[nuc] * unumpy.exp(-self.lams[nuc] * t)
            counts.append(run_count)

        print(f'Number of nuclides: {len(use_nucs)}')
        end = time.time()
        print(f'Took {round(end-start, 3)}s for {len(self.times)} steps to generate counts')
        return self.times, counts
    
    def from_groups(self, yields: list, lams: list, fissions: float,
                    fission_rate: float=None):
        counts = list()

        if type(fission_rate) is type(None):
            leading_term = fissions
            mult_term = [yields[i] * lams[i] for i in range(len(yields))]
        else:
            leading_term = fission_rate
            mult_term = yields
        for t in self.times:
            run_count = 0
            for group in range(len(yields)):
                run_count += mult_term[group] * np.exp(-lams[group] * t)
            counts.append(leading_term * run_count)
        print(f'Group delnu (summed yields): {round(sum(yields), 4)}')
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
            use_counts = np.asarray([unumpy.nominal_values(x) for x in counts[i]])
            d_counts = np.asarray([unumpy.std_devs(x) for x in counts[i]])
            plt.plot(times[i], use_counts, label=count_names[i])
            plt.fill_between(times[i], use_counts+d_counts, use_counts-d_counts, alpha=0.6)
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
                use_counts = np.asarray([unumpy.nominal_values(counts[i][j]/fissions[i]) for j in range(len(counts[i]))])
                d_counts = np.asarray([unumpy.std_devs(counts[i][j]/fissions[i]) for j in range(len(counts[i]))])
                plt.plot(times[i], use_counts, label=count_names[i])
                plt.fill_between(times[i], use_counts+d_counts, use_counts-d_counts, alpha=0.6)
            plt.legend()
            plt.xlabel('Time [s]')
            plt.ylabel('Delayed Neutron Counts per Fission')
            plt.yscale('log')
            plt.tight_layout()
            plt.savefig(f'./images/counts_per_fiss.png')
            plt.close()
        
        base_times = times[0]
        base_counts = [unumpy.nominal_values(x) for x in counts[0]]
        base_name = count_names[0]
        for i in range(len(count_names)):
            if np.any(base_times != times[i]):
                continue
            if i == 0:
                continue
            use_counts = [unumpy.nominal_values(x) for x in counts[i]]
            pcnt_diff = ((np.asarray(use_counts) - np.asarray(base_counts)) / np.asarray(base_counts) * 100)
            plt.plot(times[i], pcnt_diff, label=f'{count_names[i]}-{base_name}')
        plt.legend()
        plt.xlabel('Time [s]')
        plt.ylabel('Neutron Count Difference [\%]')
        plt.tight_layout()
        plt.savefig(f'./images/countscompare.png')
        plt.close()

        
        return

    def run_counter(self, name: str, time_list: list, count_list: list,
                   name_list: list, fission_list: list, method: str,
                   **kwargs):
        fissions = kwargs['fissions']
        try:
            fission_rate = kwargs['fission_rate']
        except KeyError:
            fission_rate = None
        if method == 'conc':
            csv_path = kwargs['csv_path']
            cutoff = kwargs['cutoff']
            times, delnu_counts = Count.from_concs(csv_path, cutoff_scale=cutoff)
        elif method == 'group':
            yields = kwargs['yields']
            lams = kwargs['lams']
            times, delnu_counts = Count.from_groups(yields, lams, fissions,
                                                    fission_rate)
        if type(fission_rate) != type(None):
            delnu = Count.get_delnu(fission_rate, delnu_counts, times)
            fission_list.append(fission_rate)
        else:
            delnu = Count.get_delnu(fissions, delnu_counts, times)
            fission_list.append(fissions)

        print(f'(Only valid if pulse) (integration of counts) {name=} {delnu=:.3E}\n')
        time_list.append(times)
        count_list.append(delnu_counts)
        name_list.append(name)
        return time_list, count_list, name_list, fission_list




if __name__ == "__main__":
    import ui
    dt = 1e-1
    tf = 500
    t0 = 0
    cutoff = 1


    time_list = list()
    count_list = list()
    name_list = list()
    fission_list = list()
    Count = DelayedCounts(dt, tf, t0)

    name = 'Static'
    csv_path = f'./results/{name}/concs.csv'
    avg_fiss_rate = 4.290E+13
    net_fiss = 1.802E+16
    fissions = avg_fiss_rate
    time_list, count_list, name_list, fission_list = Count.run_counter(name,
                                                                      time_list,
                                                                      count_list,
                                                                      name_list, 
                                                                      fission_list,
                                                                      method='conc',
                                                                      csv_path=csv_path,
                                                                      fissions=fissions,
                                                                      cutoff=cutoff)
#
#    name = 'Flowing'
#    csv_path = f'./results/{name}/concs.csv'
#    avg_fiss_rate=4.222E+13
#    net_fiss=1.773E+16
#    fissions = avg_fiss_rate
#    time_list, count_list, name_list, fission_list = Count.run_counter(name,
#                                                                      time_list,
#                                                                      count_list,
#                                                                      name_list, 
#                                                                      fission_list,
#                                                                      method='conc',
#                                                                      csv_path=csv_path,
#                                                                      fissions=fissions,
#                                                                      cutoff=cutoff)
#
#    name = 'ExFlowing'
#    csv_path = f'./results/{name}/concs.csv'
#    avg_fiss_rate=8.535e13
#    net_fiss=4.272E+16
#    fissions = avg_fiss_rate
#    time_list, count_list, name_list, fission_list = Count.run_counter(name,
#                                                                      time_list,
#                                                                      count_list,
#                                                                      name_list, 
#                                                                      fission_list,
#                                                                      method='conc',
#                                                                      csv_path=csv_path,
#                                                                      fissions=fissions,
#                                                                      cutoff=cutoff)
#    



#    name = 'Keepin Fit (Saturation)'
#    yields = [0.00063, 0.00351, 0.00310, 0.00672, 0.00211, 0.00043]
#    hls = [54.51, 21.84, 6.00, 2.23, 0.496, 0.179]
#    lams = [np.log(2)/hl for hl in hls]
#    fissions = 3E13
#    fission_rate = 3E13
#    time_list, count_list, name_list, fission_list = Count.run_counter(name,
#                                                                      time_list,
#                                                                      count_list,
#                                                                      name_list, 
#                                                                      fission_list,
#                                                                      method='group',
#                                                                      fissions=fissions,
#                                                                      cutoff=cutoff,
#                                                                      fission_rate=fission_rate,
#                                                                      yields=yields,
#                                                                      lams=lams)

#    name = 'Pulse'
#    csv_path = f'./results/{name}/concs.csv'
#    avg_fiss_rate=3.401E+17
#    net_fiss=8.502E+14
#    fissions = net_fiss
#    time_list, count_list, name_list, fission_list = Count.run_counter(name,
#                                                                      time_list,
#                                                                      count_list,
#                                                                      name_list, 
#                                                                      fission_list,
#                                                                      method='conc',
#                                                                      csv_path=csv_path,
#                                                                      fissions=fissions,
#                                                                      cutoff=cutoff)

#    name = 'Keepin Fit (Pulse)'
#    yields = [0.00063, 0.00351, 0.00310, 0.00672, 0.00211, 0.00043]
#    hls = [54.51, 21.84, 6.00, 2.23, 0.496, 0.179]
#    lams = [np.log(2)/hl for hl in hls]
#    fissions = 1E16
#    time_list, count_list, name_list, fission_list = Count.run_counter(name,
#                                                                      time_list,
#                                                                      count_list,
#                                                                      name_list, 
#                                                                      fission_list,
#                                                                      method='group',
#                                                                      csv_path=csv_path,
#                                                                      fissions=fissions,
#                                                                      cutoff=cutoff,
#                                                                      yields=yields,
#                                                                      lams=lams)

#    name = 'Static Fit'
#    # Fast Pu239
#    #yields = [0.00017, 0.00191, 0.0004, 0.00171, 0.00245, 0.00075]
#    #hls = [55.62549, 24.5045, 12.71955, 4.64925, 1.86579, 0.34046]
#    # Thermal Pu239
#    #yields = [0.00017, 0.00184, 0.0004, 0.00171, 0.00245, 0.00075]
#    #hls = [55.62548, 24.50371, 12.71953, 4.64924, 1.86579, 0.34046]
#    # Fast U235
#    yields = [0.00053, 0.00254, 0.00102, 0.0041, 0.00532, 0.00535]
#    hls = [55.6093, 24.23124, 12.71509, 4.04672, 1.64556, 0.27044]
#    # Thermal U235
#    #yields = [0.00052, 0.00238, 0.0011, 0.00386, 0.00562, 0.00541]
#    #hls = [55.63849, 24.44051, 14.40527, 4.29679, 1.71689, 0.27369]
#    lams = [np.log(2)/hl for hl in hls]
#    fissions = 1E16
#    time_list, count_list, name_list, fission_list = Count.run_counter(name,
#                                                                      time_list,
#                                                                      count_list,
#                                                                      name_list, 
#                                                                      fission_list,
#                                                                      method='group',
#                                                                      fissions=fissions,
#                                                                      cutoff=cutoff,
#                                                                      yields=yields,
#                                                                      lams=lams)
    

#    name = 'Repr Static Fit'
#    # Fast Pu239
#    # Thermal Pu239
#    # Fast U235
#    yields = [0.00047, 0.00244, 0.00102, 0.0041, 0.00532, 0.00535]
#    hls = [55.6085, 24.26064, 12.71509, 4.04672, 1.64556, 0.27044]
#    # Thermal U235
#    lams = [np.log(2)/hl for hl in hls]
#    fissions = 1E16
#    time_list, count_list, name_list, fission_list = Count.run_counter(name,
#                                                                      time_list,
#                                                                      count_list,
#                                                                      name_list, 
#                                                                      fission_list,
#                                                                      method='group',
#                                                                      fissions=fissions,
#                                                                      cutoff=cutoff,
#                                                                      yields=yields,
#                                                                      lams=lams)

#    name = 'Flowing Repr Fit'
#    # Fast U235
#    yields = [0.00045, 0.00196, 0.0011, 0.00386, 0.00562, 0.00541]
#    hls = [55.63962, 24.58418, 14.40527, 4.29679, 1.71689, 0.27369]
#    # Thermal U235
#    #yields = [0.0005, 0.00204, 0.0011, 0.00386, 0.00562, 0.00541]
#    #hls = [55.63979, 24.56228, 14.40527, 4.29679, 1.71689, 0.27369]
#    lams = [np.log(2)/hl for hl in hls]
#    fissions = 1e16
#    time_list, count_list, name_list, fission_list = Count.run_counter(name,
#                                                                      time_list,
#                                                                      count_list,
#                                                                      name_list, 
#                                                                      fission_list,
#                                                                      method='group',
#                                                                      fissions=fissions,
#                                                                      cutoff=cutoff,
#                                                                      yields=yields,
#                                                                      lams=lams)

#    name = 'Flowing Fit'
#    # Fast Pu239
#    #yields = [0.00017, 0.00183, 0.0004, 0.00171, 0.00245, 0.00075]
#    #hls = [55.31714, 23.84857, 12.71955, 4.64925, 1.86579, 0.34046]
#    # Fast U235
#    yields = [0.00051, 0.00236, 0.00102, 0.0041, 0.00532, 0.00535]
#    hls = [55.57032, 23.913, 12.71509, 4.04672, 1.64556, 0.27044]
#    # Thermal U235
#    #yields = [0.00051, 0.00207, 0.0011, 0.00386, 0.00562, 0.00541]
#    #hls = [55.63997, 24.58458, 14.40527, 4.29679, 1.71689, 0.27369]
#    lams = [np.log(2)/hl for hl in hls]
#    fissions = 1E16
#    time_list, count_list, name_list, fission_list = Count.run_counter(name,
#                                                                      time_list,
#                                                                      count_list,
#                                                                      name_list, 
#                                                                      fission_list,
#                                                                      method='group',
#                                                                      fissions=fissions,
#                                                                      cutoff=cutoff,
#                                                                      yields=yields,
#                                                                      lams=lams)

#    name = 'ExFlowing Fit'
#    # Thermal Pu239
#    yields = [0.00016, 0.0016, 0.0004, 0.00171, 0.00245, 0.00075]
#    hls = [55.63946, 24.58579, 12.71953, 4.64924, 1.86579, 0.34046]
#    lams = [np.log(2)/hl for hl in hls]
#    fissions = 1E16
#    time_list, count_list, name_list, fission_list = Count.run_counter(name,
#                                                                      time_list,
#                                                                      count_list,
#                                                                      name_list, 
#                                                                      fission_list,
#                                                                      method='group',
#                                                                      fissions=fissions,
#                                                                      cutoff=cutoff,
#                                                                      yields=yields,
#                                                                      lams=lams)

    


    Count.count_compare(count_names=name_list,
                        times=time_list,
                        counts=count_list,
                        fissions=fission_list)

