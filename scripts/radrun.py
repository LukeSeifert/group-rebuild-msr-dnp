from simple import IrradSimple
import copy
import matplotlib.pyplot as plt
import os
import numpy as np
import scipy as sp
import pandas as pd

class Run:
    """
    This class handles running the various irradiation simulations.

    """

    def __init__(self, nuc_list: list, run_omc: bool = False):
        """

        Parameters
        ----------
        nuc_list : list
            List of nuclides to analyze
        run_omc : bool, optional
            Run OpenMC irradiation simulation, by default False
        """
        self.run_omc = run_omc
        self.nuc_list = nuc_list
        self.metadict = dict()
        self.metadict['conc'] = dict()
        self.metadict['fiss'] = dict()
        self.metadict['times'] = dict()
        self.metadict['delnu'] = dict()

        simple_irrad_obj = IrradSimple(None)
        simple_irrad_obj._get_data()
        self.iaea_nucs = simple_irrad_obj.iaea_nucs
        return
    
    def _reset_metadict(self):
        self.metadict['conc'] = dict()
        self.metadict['fiss'] = dict()
        self.metadict['times'] = dict()
        self.metadict['delnu'] = dict()
        return
    
    def _simple_run(self, simple_irrad_obj: IrradSimple):
        if self.run_omc:
            simple_irrad_obj.irradiate()
        concs, times = simple_irrad_obj.collect_concentrations()
        fisses = simple_irrad_obj.collect_fissions()
        delnus = simple_irrad_obj.collect_delnu()
        version = simple_irrad_obj.name
        self.metadict['conc'][version] = concs
        self.metadict['fiss'][version] = fisses
        self.metadict['times'][version] = times
        self.metadict['delnu'][version] = delnus
        return
    
    def _fission_analysis(self):
        self.avg_fiss_rate = dict()
        for vi, version in enumerate(self.metadict['fiss'].keys()):
            fission_rates = self.metadict['fiss'][version]['net']
            times = self.metadict['times'][version]
            fissions = np.sum(fission_rates[:-1] * np.diff(times))
            target_nuc = self.irrad_objs[vi].sample_nuc
            target_fiss_rate = self.metadict['fiss'][version][target_nuc]
            target_fiss = np.sum(target_fiss_rate[:-1] * np.diff(times))
            fiss_frac_min = np.min(target_fiss / fissions)
            net_fiss = np.sum(fissions)
            avg_fiss_rate = net_fiss / times[-1]
            self.avg_fiss_rate[version] = avg_fiss_rate
            print(version)
            print(f'     {avg_fiss_rate = :.3E}')
            print(f'     {net_fiss = :.3E}')
            print(f'     {fiss_frac_min = :.3E}')
        return avg_fiss_rate, net_fiss
    
    def _delnu_analysis(self):
        for vi, version in enumerate(self.metadict['fiss'].keys()):
            delnu_data = self.metadict['delnu'][version]['net']
            final_nonzero = [i for i in reversed(delnu_data) if i > 0][0]
            delnu = final_nonzero / self.avg_fiss_rate[version]
            print(version)
            print(f'     {delnu=:.3E} from OpenMC tally')
        return

    
    def _find_max_nuc_diff(self):
        """
        Find the largest nuclide concentration difference at the final timestep
        (max L2 norm between all from first `IrradSimple` object).

        Returns
        -------
        max_norm_nuc : str
            Name of nuclide which holds largest norm
        max_norm : float
            Maximum norm of final concentration vector for each version 
        """
        versions = self.metadict['conc'].keys()
        nucs = list()
        concs = list()
        diffs = list()
        for version in versions:
            nucs += list(self.metadict['conc'][version].keys())
        nucs = list(set(nucs))
        top_nucs = 10

        for nuc in nucs:
            conc_vector = list()
            for version in versions:
                conc = self.metadict['conc'][version][nuc][-1]
                conc_vector.append(conc)
            conc_vector = np.asarray(conc_vector)
            first_element = conc_vector[0]
            percent_differences = 0.0
            if first_element > 0.0:
                percent_differences = 100 * (conc_vector - first_element) / first_element
            concs.append(conc_vector)
            conc_norm = np.linalg.norm(conc_vector)
            diffs.append(conc_norm)

        zipped_lists = list(zip(nucs, concs, diffs))
        sorted_zipped_lists = sorted(zipped_lists, key=lambda x: x[2], reverse=True)
        nucs_sorted, concs_sorted, diffs_sorted = zip(*sorted_zipped_lists)
        nucs_sorted = list(nucs_sorted)
        concs_sorted = list(concs_sorted)
        diffs_sorted = list(diffs_sorted)

        print(f'Top {top_nucs} concentration % diffs')
        for i in range(top_nucs):
            first = concs_sorted[i][0]
            pcnt_diff_vec = [100 * (i - first) / first for i in concs_sorted[i]]
            formatted_diffs = [f"{abs(diff):.3E}" for diff in pcnt_diff_vec[1:]]
            print(f'     {nucs_sorted[i]}: {formatted_diffs}')

        print(f'Top {top_nucs} DNP concentration % diffs')
        j = 0
        while j < top_nucs:
            i += 1
            if nucs_sorted[i] in self.iaea_nucs:
                first = concs_sorted[i][0]
                pcnt_diff_vec = [100 * (i - first) / first for i in concs_sorted[i]]
                formatted_diffs = [f"{abs(diff):.3E}" for diff in pcnt_diff_vec[1:]]
                print(f'     {nucs_sorted[i]}: {formatted_diffs}')
                j += 1

        return sorted_zipped_lists
    
    
    def simple_compare(self, *args):
        """
        Takes in any number of `IrradSimple` objects to be compared

        """
        self.irrad_objs = list()
        for simple_irrad_obj in args:
            self._simple_run(simple_irrad_obj) 
            self.irrad_objs.append(simple_irrad_obj)

        avgF, netF = self._fission_analysis()
        self._delnu_analysis()
        self._nuc_compare()
        self._find_max_nuc_diff()


        return avgF, netF
    

    def _nuc_compare(self):
        for nuc in self.nuc_list:
            for version in self.metadict['conc'].keys():
                conc = self.metadict['conc'][version][nuc]
                times = self.metadict['times'][version]
                plt.plot(times, conc, label=version)
            plt.xlabel('Time [s]')
            plt.ylabel('Concentration [atoms]')
            plt.legend()
            plt.tight_layout()
            try:
                plt.savefig(f'./images/{nuc}_conc.png')
            except FileNotFoundError:
                os.mkdir('./images')
                plt.savefig(f'./images/{nuc}_conc.png')
            plt.close()
        return
    
    def write_concentrations(self, simple_irrad_obj: IrradSimple):
        """
        Write the concentrations from the irradiation to a csv file.
        Number of nuclides pre-pruning: 3820
        Number after removing 0's: 1524
        Number after requiring nuc be in iaea nuc list: 309

        Parameters
        ----------
        simple_irrad_obj : IrradSimple
            Simple irradiation object
        """
        # Only non-zero and in iaea nucs
        version = simple_irrad_obj.name
        outpath = simple_irrad_obj.output_path
        self._simple_run(simple_irrad_obj)
        concs = self.metadict['conc'][version]
        final_values = {key: values[-1] for key, values in concs.items() if values[-1] > 0.0 and key in self.iaea_nucs}
        df = pd.DataFrame.from_dict(final_values, orient='index')
        output_csv_file = f'{outpath}/concs.csv'
        df.to_csv(output_csv_file, header=False)
        return


if __name__ == "__main__":
    import ui
    runner = Run(ui.nuc_list,
                 run_omc=True)
    flowing = IrradSimple(data_dict=ui.flow_data)
    static = IrradSimple(data_dict=ui.static_data)
    pulse = IrradSimple(data_dict=ui.pulse_data)
    exflow = IrradSimple(data_dict=ui.mostly_excore_data)
    reprflow = IrradSimple(data_dict=ui.flow_repr_data)
    exreprflow = IrradSimple(data_dict=ui.exflow_repr_data)
    reprstatic = IrradSimple(data_dict=ui.static_repr_data)
    #runner.simple_compare(flowing,
    #                      static,
    #                      pulse,
    #                      exflow)
    #runner.simple_compare(exflow, flowing, static, reprflow)
    #runner.write_concentrations(flowing)
    #runner.write_concentrations(static)
    #runner.write_concentrations(pulse)
    #runner.write_concentrations(exflow)
    #runner.write_concentrations(exreprflow)
    #runner.write_concentrations(reprflow)
    runner.write_concentrations(reprstatic)
