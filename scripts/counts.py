import numpy as np
from simple import IrradSimple


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
            beta_prob = self.iaea_data['  beta- %']
            Pn1 = self.iaea_data[' pn1 % ']
            Pn2 = self.iaea_data[' pn2 % ']
            Pn3 = self.iaea_data[' pn3 % ']
            prob_neutron = beta_prob * (1*Pn1 + 2*Pn2 + 3*Pn3)
            lam = self.iaea_data[' T1/2 [s] ']
            self.pns[nuc] = prob_neutron
            self.lams[nuc] = np.log(2) / lam 
        return
    
    def from_concs(self):
        self._order_iaea_data()
        
        return



if __name__ == "__main__":
    # Read in CSV data
    # Read in IAEA data
    # OR
    # Read in from group data
    # Combine
    Count = DelayedCounts(1e-3, 10)
    Count.from_concs()
