from simple import IrradSimple
import copy
import matplotlib.pyplot as plt
import os
import numpy as np
import scipy as sp

class Run:
    """
    This class handles running the various irradiation simulations.

    """

    def __init__(self, nuc_list: list, run_omc: bool = False, analyze: bool = False):
        """

        Parameters
        ----------
        nuc_list : list
            List of nuclides to analyze
        run_omc : bool, optional
            Run OpenMC irradiation simulation, by default False
        analyze : bool, optional
            Run analysis of results, by default False
        """
        self.run_omc = run_omc
        self.analyze = analyze
        self.nuc_list = nuc_list
        self.metadict = dict()
        self.metadict['conc'] = dict()
        self.metadict['fiss'] = dict()
        self.metadict['times'] = dict()
        return
    
    def _simple_run(self, simple_irrad_obj):
        if self.run_omc:
            simple_irrad_obj.irradiate()
        concs, times = simple_irrad_obj.collect_concentrations()
        fisses = simple_irrad_obj.collect_fissions()
        version = simple_irrad_obj.name
        self.metadict['conc'][version] = concs
        self.metadict['fiss'][version] = fisses
        self.metadict['times'][version] = times
        return
    
    def _fission_analysis(self):
        for vi, version in enumerate(self.metadict['fiss'].keys()):
            fission_rates = self.metadict['fiss'][version]['net']
            times = self.metadict['times'][version]
            fissions = sp.integrate.simpson(fission_rates, x=times)
            target_nuc = self.irrad_objs[vi].sample_nuc
            target_fiss_rate = self.metadict['fiss'][version][target_nuc]
            target_fiss = sp.integrate.simpson(target_fiss_rate, x=times)
            fiss_frac = target_fiss / fissions
            net_fiss = np.sum(fissions)
            avg_fiss_rate = net_fiss / times[-1]
            print(version)
            print(f'     {avg_fiss_rate=:.3E}')
            print(f'     {net_fiss=:.3E}')
        return

    
    def simple_compare(self, *args):
        """
        Takes in any number of `IrradSimple` objects to be compared

        """
        self.irrad_objs = list()
        for simple_irrad_obj in args:
            self._simple_run(simple_irrad_obj) 
            self.irrad_objs.append(simple_irrad_obj)

        if self.analyze:
            self._fission_analysis()
            self._nuc_compare()


        return
    

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





if __name__ == "__main__":
    import ui
    runner = Run(ui.nuc_list,
                 run_omc=True,
                 analyze=True)
    flowing = IrradSimple(data_dict=ui.base_case_data)
    static = IrradSimple(data_dict=ui.static_data)
    runner.simple_compare(flowing,
                          static)
