from simple import IrradSimple
import copy

class Run:
    """
    This class handles running the various irradiation simulations.

    """

    def __init__(self, run_omc=False, analyze=False):
        self.run_omc = run_omc
        self.analyze = analyze
        self.metadata = dict()
        self.metadata['conc'] = dict()
        self.metadata['fiss'] = dict()
        return
    
    def _simple_run(self, simple_irrad_obj, version):
        if self.run_omc():
            simple_irrad_obj.irradiate()
        concs = simple_irrad_obj.collect_concentrations()
        fisses = simple_irrad_obj.collect_fissions()
        self.metadata['conc'][version] = concs
        self.metadata['fiss'][version] = fisses
        return

    
    def simple_compare(self, *args):
        """
        Takes in any number of `IrradSimple` objects to be compared

        """
        for simple_irrad_obj in args:
            self._simple_run(simple_irrad_obj, version='flowing')
        return




if __name__ == "__main__":
    import ui
