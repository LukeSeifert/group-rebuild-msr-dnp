from simple import IrradSimple
import copy
import matplotlib.pyplot as plt
import os

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

    
    def simple_compare(self, *args):
        """
        Takes in any number of `IrradSimple` objects to be compared

        """
        for simple_irrad_obj in args:
            self._simple_run(simple_irrad_obj) 

        if self.analyze:
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
                 run_omc=False,
                 analyze=True)
    flowing = IrradSimple(data_dict=ui.base_case_data)
    static = IrradSimple(data_dict=ui.static_data)
    runner.simple_compare(flowing,
                          static)
