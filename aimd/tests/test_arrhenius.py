from __future__ import unicode_literals, print_function, division

"""
Created on Feb.4, 2018
"""

import unittest
import os
from aimd.diffusion import ArreheniusAnalyzer, get_conversion_factor
from pymatgen.io.vasp import Poscar

test_dir = os.path.join(os.path.abspath(os.path.dirname(__file__)),
                        'test_files')


class ArreheniusAnalyzerTest(unittest.TestCase):
    def test_arrhenius(self):
        poscar = Poscar.from_file(os.path.join(test_dir, 'arrhenius/POSCAR'))
        struct = poscar.structure
        aa = ArreheniusAnalyzer.from_csv(os.path.join(test_dir, 'arrhenius/D_T.csv'))
        temperature = 300
        specie_string = 'Li+'
        self.assertAlmostEqual(get_conversion_factor(struct, specie_string, temperature),
                               158701249.73192352, places=5)
        self.assertAlmostEqual(aa.Ea, 0.258128283041, places=5)
        self.assertAlmostEqual(aa.Ea_error, 0.0168320198494, places=5)
        self.assertAlmostEqual(aa.predict_conductivity(temperature, struct, specie_string)[0],
                               1.08036514392, places=5)
        self.assertAlmostEqual(aa.predict_conductivity(temperature, struct, specie_string)[1][0],
                               0.54509683117599272, places=5)
        self.assertAlmostEqual(aa.predict_conductivity(temperature, struct, specie_string)[1][1],
                               2.1412504667878154, places=5)
        self.assertAlmostEqual(aa.predict_diffusivity(temperature)[0],
                               6.80754024147e-09,places=5)


if __name__ == "__main__":
    unittest.main()
