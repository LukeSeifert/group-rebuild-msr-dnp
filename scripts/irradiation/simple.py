import openmc
import openmc.deplete
import numpy as np

class IrradSimple:
    """
    This class handles a simplified form of irradiation of a sample circulating
      in an MSR. This is handled by having two time based inputs, in-core and
      ex-core residence times. The sample is irradiated for the in-core time, 
      then decays over the ex-core time. This is repeated to simulate saturation
      irradiation. 
    """

    def __init__(self, t_incore_s: float, t_excore_s: float, n_MeV: float,
                 S_rate_per_s: float, batches: int,
                 output_path: str, nps: float, photon_bool: bool,
                 run_mode: str, temperature_K: float, fissile_nuc: str,
                 dens_g_cc: float, chain: str):
        """

        Parameters
        ----------
        t_incore_s : float
            The amount of time the fuel sample spends in core
        t_excore_s : float
            The amount of time the fuel sample spends out of the core
        n_MeV : float
            Neutron irradiation energy
        S_rate_per_s : float
            Irradiation source rate
        batches : int
            Number of batches to use in OpenMC
        output_path : str
            Path to write OpenMC output
        nps : float
            Number of particles to use in OpenMC simulation each batch
        photon_bool : bool
            If True, include photons
        run_mode : str
            OpenMC run mode (fixed source, plot are the relevant ones here)
        temperature_K : float
            Temperature of sample
        fissile_nuc : str
            Fissile nuclide of which sample is composed
        dens_g_cc : float
            Density of the sample in grams per cubic centimeter
        chain : str
            Path to chain xml file
        """
        self.t_incore = t_incore_s
        self.t_excore = t_excore_s
        self.n_energy = n_MeV
        self.S_rate = S_rate_per_s
        self.batches = batches
        self.output_path = output_path
        self.nps = nps
        self.photons = photon_bool
        self.run_mode = run_mode #eigenvalue', 'fixed source', 'plot', 'volume', 'particle restart'
        self.temperature = temperature_K
        self.sample_nuc = fissile_nuc
        self.sample_dens = dens_g_cc
        self.chain = chain

        self.r_outer = 10
        self.vol = 3
        self.seed = 1
        self.net_irrad_time_s = 5 * 60
        return
    
    def _settings(self):
        """
        Build the settings for OpenMC

        Returns
        -------
        settings : :class:`openmc.Settings`
        """
        settings = openmc.Settings()
        settings.batches = self.batches
        settings.output = {'path': self.output_path,
                           'summary': False,
                           'tallies': False}
        settings.particles = self.nps
        settings.photon_transport = self.photons
        settings.run_mode = self.run_mode
        settings.seed = self.seed

        source = openmc.IndependentSource()
        source.energy = openmc.stats.Uniform(self.n_energy, self.n_energy)
        source.space = openmc.stats.spherical_uniform(r_outer=self.r_outer)
        source.angle = openmc.stats.Isotropic()

        settings.source = source
        settings.temperature = {'default': self.temperature,
                                'method': 'interpolation'}
        return settings
    
    def _materials(self):
        """
        Build the materials for OpenMC

        Returns
        -------
        materials : :class:`openmc.Materials`
        """
        sample = openmc.Material()# Define the geometry
        sample.name = 'sample'
        sample.add_nuclide(self.sample_nuc, 100)
        sample.set_density('g/cc', self.sample_dens)
        sample.depletable = True
        sample.volume = self.vol

        materials = openmc.Materials([sample])
        return materials
    
    def _geometry(self):
        """
        Build the geometry for OpenMC

        Returns
        -------
        geometry : :class:`openmc.Geometry` 

        """
        sphere = openmc.Sphere(r=self.r_outer,
                               boundary_type='vacuum')
        sample_mat = self.materials[0]
        sphere_cell = openmc.Cell(region=-sphere,
                                  fill=sample_mat)
        universe = openmc.Universe(cells=[sphere_cell])
        geometry = openmc.Geometry(universe)
        return geometry
    
    def _tallies(self):
        tallies = None
        return tallies
    
    def build_model(self, xml_export=False):


        """
        Builds the OpenMC model

        Parameters
        ----------
        xml_export : bool
            If True, generate OpenMC xml files
        
        Returns
        -------
        model : :class:`openmc.model.Model`
            The built model
        
        """
        self.settings = self._settings()
        self.materials = self._materials()
        self.geometry = self._geometry()
        self.tallies = self._tallies()
        model = openmc.model.Model(self.geometry,
                                   self.materials,
                                   self.settings)#,
                                   #self.tallies)
        if xml_export:
            model.export_to_xml()
        self.model = model
        return model
        
    def _get_times_rates(self):
        """
        Get the timestep list and source rate list

        """
        self.net_irrad_time_s
        self.t_incore
        self.t_excore
        cur_t = 0
        timesteps = []
        source_rates = []

        def _step_helper(t, s, cur_t, net_t):
            break_condition = False
            timesteps.append(t)
            source_rates.append(s)
            cur_t += t
            if cur_t > net_t:
                break_condition = True
            return break_condition, cur_t


        while True:
            break_condition, cur_t = _step_helper(self.t_incore, self.S_rate, cur_t,
                                           self.net_irrad_time_s)
            if break_condition:
                break

            break_condition, cur_t = _step_helper(self.t_excore, 0, cur_t,
                                           self.net_irrad_time_s)
            if break_condition:
                break

        return timesteps, source_rates
    
    def irradiate(self):
        """
        Irradiate the sample

        """
        coupled_operator = openmc.deplete.CoupledOperator(self.model,
                                                          chain_file=self.chain,
                                                          normalization_mode='source-rate')
        timesteps, source_rates = self._get_times_rates()
        integrator = openmc.deplete.PredictorIntegrator(coupled_operator,
                                                        timesteps=timesteps,
                                                        source_rates=source_rates,
                                                        timestep_units='s')
        integrator.integrate()
        return



if __name__ == "__main__":
    t_incore_s = 10
    t_excore_s = 10
    n_MeV = 0.0253 * 1e-6
    S_rate_per_s = 1e10
    batches = 10
    output_path = './results'
    nps = 1000
    photon_bool = False
    run_mode = 'fixed source'
    temperature_K = 298.15
    fissile_nuc = 'U235'
    dens_g_cc = 10
    chain = '../../data/chain/chain_casl_pwr.xml'

    irrad = IrradSimple(t_incore_s=t_incore_s,
                        t_excore_s=t_excore_s,
                        n_MeV=n_MeV,
                        S_rate_per_s=S_rate_per_s,
                        batches=batches,
                        output_path=output_path,
                        nps=nps,
                        photon_bool=photon_bool,
                        run_mode=run_mode,
                        temperature_K=temperature_K,
                        fissile_nuc=fissile_nuc,
                        dens_g_cc=dens_g_cc,
                        chain=chain)
    model = irrad.build_model(xml_export=False)
    irrad.irradiate()