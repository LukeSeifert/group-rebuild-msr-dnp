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

    def __init__(self, t_incore_s: float):
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
        self.seed = 1
        self.temperature = temperature_K
        self.sample_nuc = fissile_nuc
        self.sample_dens = dens_g_cc

        self.s_r_outer = 10
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
        source.space = openmc.stats.spherical_uniform(r_outer=self.s_r_outer)
        source.angle = openmc.stats.Isotropic()

        settings.source = #TODO
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
        sample = openmc.Material()
        sample.name = 'sample'
        sample.add_nuclide(self.sample_nuc)
        sample.density = self.sample_dens
        sample.density_units = 'g/cm3'
        sample.depletable = True
        sample.volume = 3

        materials = openmc.Materials([sample], 100)
        return materials
    
    def _geometry(self):

    
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
        settings = self._settings()
        materials = self._materials()
        geometry = self._geometry()
        tallies = self._tallies()
        plots = self._plots()
        model = openmc.model.Model(geometry, materials, settings, tallies,
                                   plots)
        if xml_export:
            model.export_to_xml()
        return model












if __name__ == "__main__":
    pass