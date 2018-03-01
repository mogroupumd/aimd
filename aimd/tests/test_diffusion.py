from __future__ import unicode_literals, print_function, division

"""
Created on Feb.4, 2018
"""

import unittest
import os
from aimd.diffusion import DiffusivityAnalyzer, ErrorAnalysisFromDiffusivityAnalyzer
from pymatgen.core.periodic_table import Specie

test_dir = os.path.join(os.path.abspath(os.path.dirname(__file__)),
                        'test_files')


class DiffusivityAnalyzerTest(unittest.TestCase):
    def test_reading_vasprun_xml(self):
        site_distance = 3.2
        specie_string = 'Li+'
        specie = Specie.from_string(specie_string)
        vasprun_dirs = [os.path.join(test_dir, 'latp_md/RUN_{}/vasprun.xml.gz'.format(i)) for i in range(10, 30)]
        da = DiffusivityAnalyzer.from_files(vasprun_dirs, str(specie.element),
                                            spec_dict={'lower_bound': 0.5 * site_distance * site_distance,
                                                       'upper_bound': 0.5,
                                                       'minimum_msd_diff': 0.5 * site_distance * site_distance,
                                                       }
                                            )
        ea = ErrorAnalysisFromDiffusivityAnalyzer(da, site_distance=site_distance)
        summary_info = ea.get_summary_dict(oxidized_specie=specie_string)
        self.assertAlmostEqual(summary_info['diffusivity'], 7.60175023036e-05, places=5)
        self.assertAlmostEqual(summary_info['diffusivity_relative_standard_deviation'], 0.382165427856, places=5)
        self.assertAlmostEqual(summary_info['n_jump'], 100.48841284, places=5)
        self.assertAlmostEqual(summary_info['conversion_factor'], 7675326.58284, places=5)
        self.assertAlmostEqual(summary_info['temperature'], 1500, places=5)
        self.assertAlmostEqual(summary_info['conductivity'], 583.45915619213463, places=5)


if __name__ == "__main__":
    unittest.main()
