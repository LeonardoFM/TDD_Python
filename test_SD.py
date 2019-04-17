from unittest import TestCase,main
from math import *
from seismology_SD import Seismology_SD
import matplotlib.pyplot as plt

class Teste_SD(TestCase):
    """
    Test unit class for second derivative simple example using finite difference method
    """
    def setUp(self):
        self.SD = Seismology_SD(10.,100)
        self.SD.second_derivative_5pts()
        self.SD.second_derivative_3pts()
        self.SD.analytical_second_derivative_solution()
        # self.SD.error_3pts()
        self.SD.error_5pts()


    def test_uniform_distribution(self):
        """
        discretization rule: 10 discretization points per wavelength
        """
        self.assertEqual(self.SD.dx,10./(100.-1.))

    # def test_plot_analytic_wave(self):
    #     self.assertTrue(self.SD.plot_analytic_wave())

    def test_plot_comp(self):
        self.assertTrue(self.SD.plot_compar())


if __name__ == "__main__":
    main()
