import openmc
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
                 S_rate_per_s: float, batches: int, inactive: int,
                 output_path: str, nps: float, photon_bool: bool,
                 run_mode: str, temperature_K: float, fissile_nuc: str,
                 dens_g_cc: float):
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
        inactive : int
            Number of inactive batches
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
        """
        self.t_incore = t_incore_s
        self.t_excore = t_excore_s
        self.n_energy = n_MeV
        self.S_rate = S_rate_per_s
        self.batches = batches
        self.inactive = inactive
        self.output_path = output_path
        self.nps = nps
        self.photons = photon_bool
        self.run_mode = run_mode #eigenvalue', 'fixed source', 'plot', 'volume', 'particle restart'
        self.temperature = temperature_K
        self.sample_nuc = fissile_nuc
        self.sample_dens = dens_g_cc

        self.r_outer = 10
        self.vol = 3
        self.seed = 1
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
        settings.inactive = self.inactive
        settings.output = {'path': self.output_path,
                           'summary': False,
                           'tallies': False}
        settings.particles = self.nps
        settings.photon_transport = self.photons
        settings.run_mode = self.run_mode
        settings.seed = self.seed

        source = openmc.Source()
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
        sample.add_nuclide(self.sample_nuc)
        sample.density = self.sample_dens
        sample.density_units = 'g/cm3'
        sample.depletable = True
        sample.volume = self.vol

        materials = openmc.Materials([sample], 100)
        return materials
    
    def _geometry(self):
        """
        Build the geometry for OpenMC

        Returns
        -------
        geometry : :class:`openmc.Geometry` 

        """
        sphere = openmc.Sphere(r=self.r_outer,
                               boundary_type='reflective')
        sample_mat = self.materials[0]
        sphere_cell = openmc.Cell(region=-sphere,
                                  fill=sample_mat)
        universe = openmc.Universe(cells=[sphere_cell])
        geometry = openmc.Geometry(universe)
        return geometry
    
    def _tallies(self):
        tallies = None
        return tallies
    
    def _plots(self):
        plots = None
        return plots

    
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
        self.plots = self._plots()
        model = openmc.model.Model(self.geometry,
                                   self.materials,
                                   self.settings,
                                   self.tallies,
                                   self.plots)
        if xml_export:
            model.export_to_xml()
        return model


if __name__ == "__main__":
    pass