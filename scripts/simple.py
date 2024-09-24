import openmc
import openmc.deplete
import numpy as np
import shutil
import glob
import os
import pandas as pd
import time

class IrradSimple:
    """
    This class handles a simplified form of irradiation of a sample circulating
      in an MSR. This is handled by having two time based inputs, in-core and
      ex-core residence times. The sample is irradiated for the in-core time, 
      then decays over the ex-core time. This is repeated to simulate saturation
      irradiation. 
    """

    def __init__(self, data_dict: dict):
        """

        Parameters
        ----------
        data_dict : dict
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
            final_time : float
                Time at which irradiation stops
        """
        if type(data_dict) is not type(None):
            self.t_incore = data_dict['t_incore_s']
            self.t_excore = data_dict['t_excore_s']
            self.n_energy = data_dict['n_MeV']
            self.S_rate = data_dict['S_rate_per_s']
            self.batches = data_dict['batches']
            self.output_path = data_dict['output_path']
            self.nps = data_dict['nps']
            self.photons = data_dict['photon_bool']
            self.run_mode = data_dict['run_mode'] #eigenvalue', 'fixed source', 'plot', 'volume', 'particle restart'
            self.temperature = data_dict['temperature_K']
            self.sample_nuc = data_dict['fissile_nuc']
            self.sample_dens = data_dict['dens_g_cc']
            self.chain = data_dict['chain']
            self.name = data_dict['name']
            self.net_irrad_time_s = data_dict['final_time']
            self.repr_dict = data_dict['repr']

            # 3 gram spherical sample
            self.vol = 3 / self.sample_dens
            self.r_outer = (3 * self.vol / (4 * np.pi)) ** (1/3)
            self.seed = 1
            self.fiss_tally_name = 'fission-rate'

        _ = self._get_data()
        return
    
    def _format_nucs_iaea_to_omc(self, nuc_info):
        """
        Converts from 14BE to Be14
        """
        digits = ''.join(filter(str.isdigit, nuc_info[0]))
        letters = ''.join(filter(str.isalpha, nuc_info[0])) 
        if nuc_info[1] > 0:
            metastable = f'_m{nuc_info[1]}'
        else:
            metastable = -nuc_info[1] * 'negm'
        return letters.capitalize() + digits + metastable
    
    def _read_iaea_csv(self):
        """
        Read the iaea csv data and organize
        Data from https://www-nds.iaea.org/relnsd/delayedn/delayedn.html
        Download CSV evaluation and do not change
        """
        dataframe = pd.read_csv(f'../data/eval.csv',
                                skiprows=1, header=0)
        self.iaea_nucs = dataframe['nucid']
        metastable = dataframe[' liso']
        nuc_info = zip(self.iaea_nucs, metastable)
        self.iaea_nucs = [self._format_nucs_iaea_to_omc(i) for i in nuc_info]
        dataframe['nucid'] = self.iaea_nucs
        return dataframe
    
    def _check_pathing(self):
        """
        Checks output path to makes sure it doesn't already exist and creates
            it.
        """
        if os.path.exists(self.output_path):
            print(f'Path {self.output_path} already exists')
            print('Deleting existing results in 5 seconds unless interrupted')
            time.sleep(5)
            shutil.rmtree(f'{self.output_path}')
        os.makedirs(self.output_path)
        return
 
    def _get_data(self):
        """
        Generate data for use in other methods

        """
        xs_library = openmc.data.DataLibrary.from_xml('/home/luke/projects/cross-section-libraries/endfb-viii.0-hdf5/cross_sections.xml')
        self.xs_nuclide_list = [i['materials'][0] for i in xs_library.libraries][0:556]
        iaea_data = self._read_iaea_csv()
        return iaea_data
    
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
        sample.id = 1
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
        tallies_file = openmc.Tallies()

        mesh = openmc.RegularMesh()
        mesh.dimension = [1, 1, 1]
        mesh.lower_left = np.array([-self.r_outer, -self.r_outer, -self.r_outer])
        mesh.upper_right = np.array([self.r_outer, self.r_outer, self.r_outer])
        mesh_filter = openmc.MeshFilter(mesh)

        # Fission Rate tally
        fission_tally = openmc.Tally(name=self.fiss_tally_name)
        fission_tally.filters = [mesh_filter]
        fission_tally.scores = ['fission']
        fission_tally.multiply_density = True
        fission_tally.nuclides = self.xs_nuclide_list
        tallies_file.append(fission_tally)

        # Delayed Neutron Yield tally
        delnu_tally = openmc.Tally(name='delnuyield')
        delnu_tally.filters = [mesh_filter]
        delnu_tally.scores = ['delayed-nu-fission']
        delnu_tally.multiply_density = False
        delnu_tally.nuclides = list(set(self.xs_nuclide_list))# + self.iaea_nucs))
        tallies_file.append(delnu_tally)

        # Prompt Neutron Yield tally
        delnu_tally = openmc.Tally(name='pmtnuyield')
        delnu_tally.filters = [mesh_filter]
        delnu_tally.scores = ['prompt-nu-fission']
        delnu_tally.multiply_density = False
        delnu_tally.nuclides = list(set(self.xs_nuclide_list))# + self.iaea_nucs))
        tallies_file.append(delnu_tally)

        # Total Neutron Yield tally
        delnu_tally = openmc.Tally(name='nuyield')
        delnu_tally.filters = [mesh_filter]
        delnu_tally.scores = ['nu-fission']
        delnu_tally.multiply_density = False
        delnu_tally.nuclides = list(set(self.xs_nuclide_list))# + self.iaea_nucs))
        tallies_file.append(delnu_tally)

        return tallies_file
    
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
                                   self.settings,
                                   self.tallies)
        if xml_export:
            model.export_to_xml()
        self.model = model
        return model
    
    def _cleanup(self):
        """
        Moves h5 and xml files into the output path
        """
        xml_files = glob.glob('./*.xml')
        for file in xml_files:
            shutil.move(file, self.output_path)
        h5_files = glob.glob('./*.h5')
        for file in h5_files:
            shutil.move(file, self.output_path)
        return
        
    def _get_times_rates(self):
        """
        Get the timestep list and source rate list

        """
        cur_t = 0
        timesteps = []
        source_rates = []

        def _step_helper(t, s, cur_t, net_t):
            break_condition = False
            if t == 0:
                return break_condition, cur_t
            timesteps.append(t)
            source_rates.append(s)
            cur_t += t
            if cur_t >= net_t:
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
        self._check_pathing()
        self.build_model()
        coupled_operator = openmc.deplete.CoupledOperator(self.model,
                                                          chain_file=self.chain,
                                                          normalization_mode='source-rate')
        if self.repr_dict is not None:
            repr_obj = openmc.deplete.transfer_rates.TransferRates(coupled_operator,
                                                                   self.model)
            for element in list(self.repr_dict.keys()):
                repr_obj.set_transfer_rate('1', [element],
                                           transfer_rate=self.repr_dict[element])
        timesteps, source_rates = self._get_times_rates()
        integrator = openmc.deplete.PredictorIntegrator(coupled_operator,
                                                        timesteps=timesteps,
                                                        source_rates=source_rates,
                                                        timestep_units='s')
        if self.repr_dict is not None:
            integrator.transfer_rates = repr_obj
        start = time.time()
        integrator.integrate()
        end = time.time()
        print(f'Took {round(end-start, 3)} s')
        self._cleanup()
        return

    def collect_concentrations(self):
        """
        Collect concentration data from OpenMC depletion results

        Returns
        -------
        concs : dict 
            key : str
                Nuclide name
            value: :class:`np.ndarray`
                Array of concentrations over time
        times : :class:`np.ndarray`
            Time vector
        """
        concs = dict()
        results = openmc.deplete.Results(f'{self.output_path}/depletion_results.h5')
        self.times = results.get_times('s')
        self.dep_nuclides = list(results[0].index_nuc.keys())
        for nuc in self.dep_nuclides:
            _, nuc_conc = results.get_atoms('1',
                                            nuc=nuc)
            concs[nuc] = nuc_conc
        return concs, self.times
    
    def collect_delnu(self):
        delnu = dict()
        for i in range(len(self.times)):
            sp = openmc.StatePoint(f'{self.output_path}/openmc_simulation_n{i}.h5')
            for tally in sp.tallies.keys():
                tally_data = sp.get_tally(id=tally)
                if 'delnuyield' in tally_data.name:
                    delnu['d'] = dict()
                    use_dict = delnu['d']
                    df = tally_data.get_pandas_dataframe(filters=False, scores=False, derivative=False, paths=False)
                    try:
                        df['mean'] = df['mean'] * self.S_rate
                    except KeyError:
                        continue
                    df_sorted = df.sort_values(by='mean', ascending=False)
                    df_sorted = df_sorted.reset_index(drop=True)

                    if i == 0:
                        for nuc in df_sorted['nuclide']:
                            use_dict[nuc] = np.zeros(len(self.times))
                        use_dict['net'] = np.zeros(len(self.times))
                    for nuc_i, nuc in enumerate(df_sorted['nuclide']):
                        use_dict[nuc][i] = df_sorted['mean'][nuc_i]
                        use_dict['net'][i] += use_dict[nuc][i]

                if 'pmtnuyield' in tally_data.name:
                    delnu['p'] = dict()
                    use_dict = delnu['p']
                    df = tally_data.get_pandas_dataframe(filters=False, scores=False, derivative=False, paths=False)
                    try:
                        df['mean'] = df['mean'] * self.S_rate
                    except KeyError:
                        continue
                    df_sorted = df.sort_values(by='mean', ascending=False)
                    df_sorted = df_sorted.reset_index(drop=True)

                    if i == 0:
                        for nuc in df_sorted['nuclide']:
                            use_dict[nuc] = np.zeros(len(self.times))
                        use_dict['net'] = np.zeros(len(self.times))
                    for nuc_i, nuc in enumerate(df_sorted['nuclide']):
                        use_dict[nuc][i] = df_sorted['mean'][nuc_i]
                        use_dict['net'][i] += use_dict[nuc][i]

        fiss_keys = list(delnu.keys())
        for nuc in fiss_keys:
            if np.all(delnu[nuc] <= 1e-12 * np.ones(len(self.times))):
                del delnu[nuc]
        return delnu
    
    def collect_fissions(self, print_fiss=False):
        """
        Collect fission tally data from OpenMC depletion results

        Returns
        -------
        fissions : dict
            key : str
                Nuclide name
            value : :class:`np.ndarray`
                Array of net fissions over time
        """
        fissions = dict()
        for i in range(len(self.times)):
            sp = openmc.StatePoint(f'{self.output_path}/openmc_simulation_n{i}.h5')
            for tally in sp.tallies.keys():
                tally_data = sp.get_tally(id=tally)
                if self.fiss_tally_name in tally_data.name:
                    df = tally_data.get_pandas_dataframe(filters=False, scores=False, derivative=False, paths=False)
                    try:
                        df['mean'] = df['mean'] * self.S_rate
                    except KeyError:
                        continue
                    df_sorted = df.sort_values(by='mean', ascending=False)
                    df_sorted = df_sorted.reset_index(drop=True)
                    if i == 0:
                        for nuc in df_sorted['nuclide']:
                            fissions[nuc] = np.zeros(len(self.times))
                        fissions['net'] = np.zeros(len(self.times))
                    for nuc_i, nuc in enumerate(df_sorted['nuclide']):
                        fissions[nuc][i] = df_sorted['mean'][nuc_i]
                        fissions['net'][i] += fissions[nuc][i]

        fiss_keys = list(fissions.keys())
        for nuc in fiss_keys:
            if np.all(fissions[nuc] <= 1e-12 * np.ones(len(self.times))):
                del fissions[nuc]
        if print_fiss:
            print(fissions)
        return fissions


if __name__ == "__main__":
    import ui

    irrad = IrradSimple(data_dict=ui.pulse_data)
    #irrad.irradiate()
    irrad.collect_concentrations()
    irrad.collect_fissions()
    irrad.collect_delnu()
