import numpy as np
from math import *
import matplotlib.pyplot as plt

class Sismology_FD():

    def __init__(self,xmax,nx):
        #initial variables
        self.xmax = xmax                          #physical domain
        self.nx = nx                             #number of samples
        self.dx = self.grid_dx()                  #grid increment
        self.x = np.linspace(0,self.xmax,self.nx) #space coordinates
        # Analytical solution
        self.l = self.wavelength()
        self.k = self.wavenumber()               #wavenumber
        self.f = analytical_solution()

    def grid_dx(self):
        """
        this function create a uniform distribution
        """
        return self.xmax/(self.nx-1)

    def wavelength(self):
        """
        this function respect the number of dofs per wavelength
        """
        return 20*self.dx

    def wavenumber(self):
        """
        numerical parameter k for analytical solution
        """
        return 2*pi/self.l

    def analytical_solution(self):
        return np.sin(self.k*self.x)
