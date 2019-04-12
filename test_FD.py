from unittest import TestCase,main
from math import *
from sismology_FD import Sismology_FD

class Teste_FD(TestCase):
    """
    Test unit class for simple simology example using finite difference method
    """
    def setUp(self):
        self.FD = Sismology_FD(10.,200)

    def test_uniform_distribution(self):
        self.assertEqual(self.FD.dx,10./(200.-1.))

    def test_wavelength_imposition(self):
        self.assertEqual(self.FD.l,20*(10./(200.-1.)))

    def test_wavenumber(self):
        self.assertEqual(self.FD.k,2*pi/(20*(10./(200.-1.))))


if __name__ == "__main__":
    main()
