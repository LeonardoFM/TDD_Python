from unittest import TestCase,main
from math import *
from sismology_FD import Sismology_FD
import matplotlib.pyplot as plt

class Teste_FD(TestCase):
    """
    Test unit class for first derivative simple example using finite difference method (by simology course)
    """
    def setUp(self):
        self.FD = Sismology_FD(10.,200)
        self.FD.firt_derivative()
        self.FD.analytical_derivative_solution()
        self.FD.error()


    def test_uniform_distribution(self):
        """
        discretization rule: 10 discretization points per wavelength
        """
        self.assertEqual(self.FD.dx,10./(200.-1.))

    def test_wavelength_imposition(self):
        """
        if we have 200 points than the wavelength rule is 20*dx
        """
        self.assertEqual(self.FD.l,20*(10./(200.-1.)))

    def test_wavenumber(self):
        """
        this parameter k is the "number of thumb" for numerical approx.
        """
        self.assertEqual(self.FD.k,2*pi/(20*(10./(200.-1.))))

    # def test_plot_analytic_wave(self):
    #      self.assertTrue(self.FD.plot_analytic_wave())

    def test_plot_comp(self):
         self.assertTrue(self.FD.plot_compar())

    # def test_plot_comp(self):
    #     self.assertTrue(self.FD.plot_number_dof_per_wavelength())


if __name__ == "__main__":
    main()
