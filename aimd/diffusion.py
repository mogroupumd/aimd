# coding: utf-8
# Copyright (c) MoGroup at UMD.
# Distributed under the terms of the MIT License.

from __future__ import division, unicode_literals
import numpy as np
from monty.json import MSONable
from scipy import stats
from scipy.optimize import curve_fit
from pymatgen.io.vasp.outputs import Vasprun
from pymatgen.util.coord import pbc_diff
from pymatgen.core.structure import Structure
from pymatgen.core.periodic_table import Specie
import csv

__author__ = "Xingfeng He"
__email__ = "albert2@umd.edu"
__version__ = "0.1"
__date__ = "6/6/2017"
__credit__ = "Pymatgen Development Team"


class DiffusivityAnalyzer(MSONable):
    def __init__(self, structure, displacements, specie, temperature,
                 time_step, step_skip, time_intervals_number=1000,
                 spec_dict=None):
        """
        Calculate MSD from pre-processed data, and implemented linear fitting to obtain diffusivity.

        :param structure (Structure): initial structure
        :param displacements (np.array): numpy array, shape is [n_ions, n_steps, axis]
        :param specie (str): species string, can be Li or Li+, make sure structure has oxidation
                state accordingly.
        :param temperature (float): temperature of MD
        :param time_step (float): time step in MD
        :param step_skip (int): Sampling frequency of the displacements (
                time_step is multiplied by this number to get the real time
                between measurements)
        :param time_intervals_number (int): number of time intervals. Default is 1000
                means there are ~1000 time intervals.
        :param spec_dict (dict): spec dict of linear fitting. Default is
                {'lower_bound': 4.5, 'upper_bound': 0.5, 'minimum_msd_diff': 4.5}
                lower_bound is in unit of Angstrom square
                upper_bound is in unit of total time. 0.5 means upper fitting bound is 0.5*t_total
                minimum_msd_diff is in unit of Angstrom square. msd[upper_bound] - msd[lower_bound] should larger
                than minimum_msd_diff to do linear fitting.
        """
        spec_dict = spec_dict if spec_dict is not None else {'lower_bound': 4.5, 'upper_bound': 0.5,
                                                             'minimum_msd_diff': 4.5}
        if not {'lower_bound', 'upper_bound', 'minimum_msd_diff'} <= set(spec_dict.keys()):
            raise Exception("spec_dict does not have enough parameters.")
        time_step_displacements = time_step * step_skip
        # prepare
        indices = []
        framework_indices = []
        for i, site in enumerate(structure):
            if site.specie.symbol == specie:
                indices.append(i)
            else:
                framework_indices.append(i)
        if len(indices) == 0:
            raise Exception("There is no specie {} in the structure".format(specie))
        if len(framework_indices) == 0:
            dc = displacements
        else:
            framework_disp = displacements[framework_indices]
            drift = np.average(framework_disp, axis=0)[None, :, :]
            dc = displacements - drift
        df = structure.lattice.get_fractional_coords(dc)
        displacements_final_diffusion_ions = dc[indices]
        displacements_frac_final_diffusion_ions = df[indices]
        n_ions, n_steps, dim = displacements_final_diffusion_ions.shape

        # time intervals, dt
        dt_indices = np.arange(1, n_steps, max(int((n_steps - 1) / time_intervals_number), 1))
        dt = dt_indices * time_step_displacements

        # calculate msd
        # define functions, algorithm from
        # http://stackoverflow.com/questions/34222272/computing-mean-square-displacement-using-python-and-fft
        def autocorrelation_fft(x):
            N = x.shape[0]
            F = np.fft.fft(x, n=2 * N)
            PSD = F * F.conjugate()
            res = np.fft.ifft(PSD)
            res = (res[:N]).real
            n = N * np.ones(N) - np.arange(N)
            return res / n

        def one_ion_msd_fft(r, dt_indices):
            """
            r (np.array, shape is typically [n_step,3], n_step is number of steps, 3 is 3 dimentions)
            """
            # ------------ S1
            n_step, dim = r.shape
            r_square = np.square(r)
            r_square = np.append(r_square, np.zeros((1, dim)), axis=0)  # (n_step+1, 3)
            S1_component = np.zeros((dim, n_step))  # (dim, n_step)
            r_square_sum = 2 * np.sum(r_square, axis=0)  # (3)
            for i in range(n_step):
                r_square_sum = r_square_sum - r_square[i - 1, :] - r_square[n_step - i, :]
                S1_component[:, i] = r_square_sum / (n_step - i)
            S1 = np.sum(S1_component, axis=0)

            # ------------ S2
            S2_component = np.array([autocorrelation_fft(r[:, i]) for i in range(r.shape[1])])  # (dim, N)
            S2 = np.sum(S2_component, axis=0)

            # ------------ return
            return (S1 - 2 * S2)[dt_indices], (S1_component - 2 * S2_component)[:, dt_indices]

        n_dt = len(dt_indices)
        msd_by_ions = np.empty([0, n_dt])  # shape of n_ions * n_dt
        msd_component_by_ions = np.empty([3, 0, n_dt])  # shape of 3 * n_ions * n_dt
        for i in range(n_ions):
            msd_i, msd_component_i = one_ion_msd_fft(displacements_final_diffusion_ions[i, :, :], dt_indices)
            msd_by_ions = np.append(msd_by_ions,
                                    msd_i.reshape(1, n_dt),
                                    axis=0)
            msd_component_by_ions = np.append(msd_component_by_ions,
                                              msd_component_i.reshape(3, 1, n_dt),
                                              axis=1)
        msd = np.average(msd_by_ions, axis=0)
        msd_component = np.average(msd_component_by_ions, axis=1)

        # further things, 1. determine lower_index, upper_index 2. linear fitting, 3. error bar

        # one headache, how about error in different axis
        lower_bound_index = len(msd[msd < spec_dict['lower_bound']])
        upper_bound_index = int(len(msd) * spec_dict['upper_bound']) - 1
        """
        if lower_bound_index >= upper_bound_index - 2:
            raise Exception("Maximum MSD is {:.2f}. ".format(max(msd)) + \
                            "MSD array has shape of {}. ".format(msd.shape) + \
                            "Lower bound index is {}, upper bound index is {}. ".format(lower_bound_index,
                                                                                        upper_bound_index) + \
                            "There is no enough data to fit. " + \
                            "Please consider extending your MD simulation or increasing the temperature.")

        if msd[upper_bound_index] - msd[lower_bound_index] < spec_dict['minimum_msd_diff']:
            raise Exception(
                "Maximum MSD is {:.2f}. ".format(max(msd)) + \
                "MSD at lower bound is {:.2f}, MSD at upper bound is {:.2f}. The MSD fitting range is too small. " \
                .format(msd[lower_bound_index], msd[upper_bound_index]) + \
                "Please consider extending your MD simulation or increasing the temperature.")
        """
        if lower_bound_index >= upper_bound_index - 2 or \
                                msd[upper_bound_index] - msd[lower_bound_index] < spec_dict['minimum_msd_diff']:
            slope = -1
            intercept = -1
            slope_components = np.zeros(dim)
        else:
            slope, intercept, _, _, _ = stats.linregress(dt[lower_bound_index:upper_bound_index + 1],
                                                         msd[lower_bound_index:upper_bound_index + 1])
            slope_components = np.zeros(dim)
            for i in range(dim):
                s, _, _, _, _ = stats.linregress(dt[lower_bound_index:upper_bound_index + 1],
                                                 msd_component[i, :][lower_bound_index:upper_bound_index + 1])
                slope_components[i] = s

        self.structure = structure
        self.indices = indices
        self.framework_indices = framework_indices
        self.disp = displacements
        self.displacements_final_diffusion_ions = displacements_final_diffusion_ions
        self.specie = specie
        self.temperature = temperature
        self.time_step = time_step
        self.step_skip = step_skip
        self.time_step_displacements = time_step_displacements
        self.time_intervals_number = time_intervals_number
        self.spec_dict = spec_dict
        if len(framework_indices) == 0:
            self.max_framework_displacement = 0.0
        else:
            self.max_ion_displacements = np.max(np.sum(
                dc ** 2, axis=-1) ** 0.5, axis=1)
            self.max_framework_displacement = \
                np.max(self.max_ion_displacements[framework_indices])

        self.dt = dt
        self.lower_bound = spec_dict['lower_bound']
        self.upper_bound = spec_dict['upper_bound']
        self.lower_bound_index = lower_bound_index
        self.upper_bound_index = upper_bound_index

        self.msd = msd
        self.msd_by_ions = msd_by_ions
        self.msd_component = msd_component
        self.diffusivity = slope / (20 * dim)
        self.diffusivity_components = slope_components / 20

    def get_summary_dict(self, oxidized_specie=None):
        """
        A summary of information
        :param oxidized_specie (str): specie string with oxidation state. If provided or specie in initial
            function is oxidized, it will calculate conductivity based on nernst-einstein relationship.
        :return: dict of diffusion information
            keys: D, D_components, specie, step_skip, temperature, msd, msd_component, dt, time_intervals_number
                  spec_dict
        """
        d = {"diffusivity": self.diffusivity,
             "diffusivity_components": self.diffusivity_components,
             "specie": self.specie,
             "step_skip": self.step_skip,
             "temperature": self.temperature,
             "msd": self.msd,
             "msd_component": self.msd_component,
             "dt": self.dt,
             "time_intervals_number": self.time_intervals_number,
             "spec_dict": self.spec_dict}
        oxi = False
        if oxidized_specie:
            df_sp = Specie.from_string(oxidized_specie)
            oxi = True
        else:
            try:
                df_sp = Specie.from_string(self.specie)
                oxi = True
            except:
                pass
        if oxi:
            factor = get_conversion_factor(self.structure, df_sp, self.temperature)
            d['conductivity'] = factor * self.diffusivity
            d['conductivity_components'] = factor * self.diffusivity_components
            d['conversion_factor'] = factor
            d['oxidation_state'] = df_sp.oxi_state
        return d

    @classmethod
    def from_structures(cls, structures, specie, temperature,
                        time_step, step_skip, time_intervals_number=1000,
                        spec_dict=None):
        """
        Convenient constructor that takes in a list of Structure objects to
        perform diffusion analysis.

        :param structures ([Structure]): list of Structure objects:
        :param specie (str): species string, like Li, Li+
        :param temperature (float): temperature of MD
        :param time_step (float): time step in MD
        :param step_skip (int): Sampling frequency of the displacements (
                time_step is multiplied by this number to get the real time
                between measurements)
        :param time_intervals_number (int): number of time intervals. Default is 1000
                means there are ~1000 time intervals.
        :param spec_dict (dict): spec dict of linear fitting. Default is
                {'lower_bound': 4.5, 'upper_bound': 0.5, 'minimum_msd_diff': 4.5}
                lower_bound is in unit of Angstrom square
                upper_bound is in unit of total time. 0.5 means upper fitting bound is 0.5*t_total
                minimum_msd_diff is in unit of Angstrom square. msd[upper_bound] - msd[lower_bound] should larger
                than minimum_msd_diff to do linear fitting.
        """
        p = []
        for i, s in enumerate(structures):
            if i == 0:
                structure = s
            p.append(np.array(s.frac_coords)[:, None])

        p.insert(0, p[0])
        p = np.concatenate(p, axis=1)
        dp = p[:, 1:] - p[:, :-1]
        dp = dp - np.round(dp)
        f_disp = np.cumsum(dp, axis=1)
        disp = structure.lattice.get_cartesian_coords(f_disp)

        return cls(structure, disp, specie, temperature,
                   time_step, step_skip=step_skip,
                   time_intervals_number=time_intervals_number,
                   spec_dict=spec_dict)

    @classmethod
    def from_vaspruns(cls, vaspruns, specie,
                      time_intervals_number=1000,
                      spec_dict=None):
        """
        Convenient constructor that takes in a list of Vasprun objects to
        perform diffusion analysis.

        :param vaspruns ([Vasprun]): List of Vaspruns (ordered):
        :param specie (str): species string, like Li, Li+
        :param time_intervals_number (int): number of time intervals. Default is 1000
                means there are ~1000 time intervals.
        :param spec_dict (dict): spec dict of linear fitting. Default is
                {'lower_bound': 4.5, 'upper_bound': 0.5, 'minimum_msd_diff': 4.5}
                lower_bound is in unit of Angstrom square
                upper_bound is in unit of total time. 0.5 means upper fitting bound is 0.5*t_total
                minimum_msd_diff is in unit of Angstrom square. msd[upper_bound] - msd[lower_bound] should larger
                than minimum_msd_diff to do linear fitting.
        """

        def get_structures(vaspruns):
            for i, vr in enumerate(vaspruns):
                if i == 0:
                    step_skip = vr.ionic_step_skip or 1
                    final_structure = vr.initial_structure
                    temperature = vr.parameters['TEEND']
                    time_step = vr.parameters['POTIM']
                    yield step_skip, temperature, time_step
                # check that the runs are continuous
                fdist = pbc_diff(vr.initial_structure.frac_coords,
                                 final_structure.frac_coords)
                if np.any(fdist > 0.001):
                    raise ValueError('initial and final structures do not '
                                     'match.')
                final_structure = vr.final_structure

                assert (vr.ionic_step_skip or 1) == step_skip
                for s in vr.ionic_steps:
                    yield s['structure']

        s = get_structures(vaspruns)
        step_skip, temperature, time_step = next(s)

        return cls.from_structures(structures=s, specie=specie,
                                   temperature=temperature, time_step=time_step, step_skip=step_skip,
                                   time_intervals_number=time_intervals_number, spec_dict=spec_dict)

    @classmethod
    def from_files(cls, filepaths, specie, step_skip=10, ncores=None,
                   time_intervals_number=1000,
                   spec_dict=None):
        """
        Convenient constructor that takes in a list of vasprun.xml paths to
        perform diffusion analysis.

        :param filepaths ([str]): List of paths to vasprun.xml files of runs, ordered.
        :param specie (str): species string, like Li, Li+
        :param step_skip (int): Sampling frequency of the displacements (
                time_step is multiplied by this number to get the real time
                between measurements)
        :param ncores (int): Numbers of cores to use for multiprocessing. Can
                speed up vasprun parsing considerably. Defaults to None,
                which means serial. It should be noted that if you want to
                use multiprocessing, the number of ionic steps in all vasprun
                .xml files should be a multiple of the ionic_step_skip.
                Otherwise, inconsistent results may arise. Serial mode has no
                such restrictions.
        :param time_intervals_number (int): number of time intervals. Default is 1000
                means there are ~1000 time intervals.
        :param spec_dict (dict): spec dict of linear fitting. Default is
                {'lower_bound': 4.5, 'upper_bound': 0.5, 'minimum_msd_diff': 4.5}
                lower_bound is in unit of Angstrom square
                upper_bound is in unit of total time. 0.5 means upper fitting bound is 0.5*t_total
                minimum_msd_diff is in unit of Angstrom square. msd[upper_bound] - msd[lower_bound] should larger
                than minimum_msd_diff to do linear fitting.
        """
        if ncores is not None and len(filepaths) > 1:
            import multiprocessing
            p = multiprocessing.Pool(ncores)
            vaspruns = p.imap(_get_vasprun,
                              [(fp, step_skip) for fp in filepaths])
            analyzer = cls.from_vaspruns(vaspruns, specie=specie,
                                         time_intervals_number=time_intervals_number,
                                         spec_dict=spec_dict)
            p.close()
            p.join()
            return analyzer
        else:
            def vr(filepaths):
                offset = 0
                for p in filepaths:
                    v = Vasprun(p, ionic_step_offset=offset,
                                ionic_step_skip=step_skip)
                    yield v
                    # Recompute offset.
                    offset = (-(v.nionic_steps - offset)) % step_skip

            return cls.from_vaspruns(vr(filepaths), specie=specie,
                                     time_intervals_number=time_intervals_number,
                                     spec_dict=spec_dict)

    def as_dict(self):
        return {
            "@module": self.__class__.__module__,
            "@class": self.__class__.__name__,
            "structure": self.structure.as_dict(),
            "displacements": self.disp.tolist(),
            "specie": self.specie,
            "temperature": self.temperature,
            "time_step": self.time_step,
            "step_skip": self.step_skip,
            "time_intervals_number": self.time_intervals_number,
            "spec_dict": self.spec_dict}

    @classmethod
    def from_dict(cls, d):
        structure = Structure.from_dict(d["structure"])
        return cls(structure, np.array(d["displacements"]), specie=d["specie"],
                   temperature=d["temperature"], time_step=d["time_step"],
                   step_skip=d["step_skip"], time_intervals_number=d["time_intervals_number"],
                   spec_dict=d['spec_dict'])


class ErrorAnalysisFromDiffusivityAnalyzer(object):
    def __init__(self, diffusivity_analyzer, site_distance=3.0):
        """
        Estimate the relative standard deviation (RSD) of D from the equation:
        RSD = 3.43/sqrt(N_jump) + 0.04

        :param diffusivity_analyzer (DiffusivityAnalyzer object):
        :param site_distance (float): the site distance between diffusion ions (averagely)
        """
        n_jump = len(diffusivity_analyzer.indices) * \
                 np.max(diffusivity_analyzer.msd) / (site_distance * site_distance)
        n_jump_component = len(diffusivity_analyzer.indices) * \
                           np.max(diffusivity_analyzer.msd_component, axis=1) / (site_distance * site_distance)
        RSD_D = 3.43 / np.sqrt(n_jump) + 0.04
        RSD_D_component = [None, None, None]
        for i in range(3):
            RSD_D_component[i] = 3.43 / np.sqrt(n_jump_component[i]) + 0.04
        self.diffusivity_analyzer = diffusivity_analyzer
        self.n_jump = n_jump
        self.n_jump_component = n_jump_component
        self.RSD_D = RSD_D
        self.RSD_D_component = np.array(RSD_D_component)

    def get_summary_dict(self, oxidized_specie=None):
        """
        A summary of information
        :param oxidized_specie (str): specie string with oxidation state. If provided or specie in initial
            function is oxidized, it will calculate conductivity based on nernst-einstein relationship.
        :return: dict of diffusion information
        """
        d = self.diffusivity_analyzer.get_summary_dict(oxidized_specie=oxidized_specie)
        d['n_jump'] = self.n_jump
        d['n_jump_component'] = self.n_jump_component
        d['diffusivity_relative_standard_deviation'] = self.RSD_D
        d['diffusivity_standard_deviation'] = self.RSD_D * d['diffusivity']
        d['diffusivity_component_relative_standard_deviation'] = self.RSD_D_component
        d['diffusivity_component_relative_standard_deviation'] = self.RSD_D_component * d['diffusivity_components']
        return d


def _get_vasprun(args):
    """
    Internal method to support multiprocessing.
    """
    return Vasprun(args[0], ionic_step_skip=args[1],
                   parse_dos=False, parse_eigen=False)


class ArreheniusAnalyzer(object):
    def __init__(self, temperatures, diffusivities, diffusivity_errors=None):
        """
        Fitting arrehenius relationship from temperatures, diffusivities and diffusivity_error
        :param temperatures (List): list of temperatures
        :param diffusivities (List): list of diffusivities at different temperatures
        :param diffusivity_errors (List): optional, list of diffusivity error at different temperatures
        """
        slope_to_eV = -8.617e-5 * 1000 * np.log(10)

        def linear(x, k, b):
            return k * x + b

        x = np.array([1000.0 / i for i in temperatures])
        y = np.array([np.log10(i) for i in diffusivities])
        if diffusivity_errors is None:
            [slope, intercept], cov = curve_fit(linear, x, y)
            slope_sigma = np.sqrt(np.diag(cov))[0]
            intercept_sigma = np.sqrt(np.diag(cov))[1]
            y_error = None
        else:
            y_error = [np.log10(np.e) * diffusivity_errors[i] / diffusivities[i] for i in range(len(diffusivities))]
            [slope, intercept], cov = curve_fit(linear, x, y, sigma=y_error, absolute_sigma=True)
            slope_sigma = np.sqrt(np.diag(cov))[0]
            intercept_sigma = np.sqrt(np.diag(cov))[1]

        self.temperatures = temperatures
        self.diffusivities = diffusivities
        self.x = x  # 1000/T
        self.y = y  # log10(D)
        self.diffusivity_errors = diffusivity_errors
        self.y_error = y_error
        self.Ea = slope_to_eV * slope
        self.Ea_error = -1 * slope_to_eV * slope_sigma
        self.intercept = intercept
        self.intercept_sigma = intercept_sigma
        self.slope = slope
        self.slope_sigma = slope_sigma

    def predict_diffusivity(self, temperature):
        """

        :param temperature (float): target temperature
        :return: corresponding diffusivity, and diffusivity range based on error of linear fitting
        """
        logD = self.slope * (1000.0 / temperature) + self.intercept
        logD_sigma = np.sqrt(np.power(self.slope_sigma * (1000.0 / temperature), 2) + \
                             np.power(self.intercept_sigma, 2))
        D_min = np.power(10, logD - logD_sigma)
        D_max = np.power(10, logD + logD_sigma)
        return np.power(10, logD), [D_min, D_max]

    def predict_conductivity(self, temperature, structure, specie):
        """

        :param temperature (float): target temperature
        :param structure (Structure): one structure, used to calculate convertion factor from diffusivity
            to conductivity
        :param specie (str): string of diffusion specie, should contain oxidation state, such as Li+, O2-
        :return: corresponding conductivity, and conductivity range based on error of linear fitting
        """
        D, [D_min, D_max] = self.predict_diffusivity(temperature)
        factor = get_conversion_factor(structure, specie, temperature)
        return factor * D, [factor * D_min, factor * D_max]

    def get_arrhenius_plot(self):
        from pymatgen.util.plotting import pretty_plot
        plt = pretty_plot(12, 8)
        arr = np.power(10, self.slope * self.x + self.intercept)
        plt.plot(self.x, self.diffusivities, 'ko', self.x, arr, 'k--', markersize=10)
        plt.errorbar(self.x, self.diffusivities, yerr=self.diffusivity_errors,
                     fmt='ko', ecolor='k', capthick=2, linewidth=2)
        ax = plt.axes()
        ax.set_yscale('log')
        plt.text(0.6, 0.85, "E$_a$ = {:.3f} eV".format(self.Ea),
                 fontsize=30, transform=plt.axes().transAxes)
        plt.ylabel("D (cm$^2$/s)")
        plt.xlabel("1000/T (K$^{-1}$)")
        plt.tight_layout()
        return plt

    @classmethod
    def from_csv(cls, csv_file):
        with open(csv_file, 'r') as csvf:
            a = csv.reader(csvf, delimiter=str(","))
            data_keys = []
            data_list = []
            for row_index, row in enumerate(a):
                if row_index == 0:
                    data_keys = row
                    data_list = [[] for _ in range(len(row))]
                else:
                    for col_index, col in enumerate(row):
                        data_list[col_index].append(float(col))
        data_dict = {data_keys[i]: data_list[i] for i in range(len(data_keys))}
        if set(data_keys) >= set(['T', 'D', 'D_error']):
            return cls(data_dict['T'], data_dict['D'], data_dict['D_error'])
        elif set(data_keys) >= set(['T', 'D']):
            return cls(data_dict['T'], data_dict['D'])
        else:
            raise Exception("Please make sure the header name in csv file is T, D, D_error(optional)")


def get_conversion_factor(structure, specie, temperature):
    """
    Conversion factor to convert between cm^2/s diffusivity measurements and
    mS/cm conductivity measurements based on number of atoms of diffusing
    species.
    :param structure (Structure): Input structure.
    :param specie (string/specie): Diffusing species string, must contain oxidation state.
    :param temperature (float): Temperature of the diffusion run in Kelvin.
    :return: Conversion factor.
        Conductivity (in mS/cm) = Conversion Factor * Diffusivity (in cm^2/s)
    """
    if type(specie) is Specie:
        df_sp = specie
    else:
        try:
            df_sp = Specie.from_string(specie)
        except:
            raise Exception("Please provide oxidation decorated specie, like Li+, O2-")
    z = df_sp.oxi_state
    el, occu = structure.composition.items()[0]
    if isinstance(el, Specie):  # oxidation decorated structure
        n = structure.composition[specie]
    else:
        n = structure.composition[str(df_sp.element)]
    if n == 0:
        raise Exception("No specie {} in the structure composition: {}".format(specie, structure.composition))
    vol = structure.volume * 1e-24  # units cm^3
    N_A = 6.022140857e+23
    e = 1.6021766208e-19
    R = 8.3144598
    return 1000 * n / (vol * N_A) * z ** 2 * (N_A * e) ** 2 \
           / (R * temperature)
