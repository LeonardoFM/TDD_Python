import numpy as np
from math import *
import matplotlib.pyplot as plt

class Seismology_FD():
    """
    Example for first derivative using finite difference method (by simology course)

    https://www.coursera.org/learn/computers-waves-simulations/

    """
    def __init__(self,xmax,nx):
        #initial variables
        self.xmax = xmax                          #physical domain
        self.nx = nx                             #number of samples
        self.dx = self.grid_dx()                  #grid increment
        self.x = np.linspace(0,self.xmax,self.nx) #space coordinates
        #initiation of numerical and analytical derivatives
        self.nder = np.zeros(self.nx) #numerical derivative
        self.ader = np.zeros(self.nx) #anal
        # Analytical solution
        self.l = self.wavelength()
        self.k = self.wavenumber()               #wavenumber
        self.f = self.analytical_solution()

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
        """
        simple sine wave
        """
        return np.sin(self.k*self.x)

    def analytical_derivative_solution(self):
        """
        analytical derivative of given function
        """
        self.ader = self.k*np.cos(self.k*self.x)
        #exclude boundaries
        self.ader[0] = 0.0
        self.ader[self.nx-1] = 0.0


    def firt_derivative(self):
        """
        first derivative with two points
        """

        #numerical derivative of a given Function
        for i in range(1,self.nx-1):
            self.nder[i] = (self.f[i+1] - self.f[i-1])/(2*self.dx)

    def error(self):
        """
        class error calculation: max norm
        """
        self.rms = np.sqrt(np.mean(self.nder-self.ader)**2)

    def plot_analytic_wave(self):
        try:
            plt.plot(self.x,self.f)
            plt.title("sin function")
            plt.xlabel("x, m")
            plt.ylabel("Amplitude")
            plt.xlim((0,self.xmax))
            plt.grid()
            plt.show()
            return True
        except:
            print("An exception occurred ")
        else:
            return False

    def plot_compar(self):
        try:
            plt.plot(self.x,self.nder,label="Numerical derivative, 2 points",marker="+",color="blue")
            plt.plot(self.x,self.ader,label="Analytical derivative",lw=2,ls="-",color="black")
            plt.plot(self.x,self.nder-self.ader,label="Diference",lw=2,ls=":")
            plt.title(f"First derivative, error = {self.rms}")
            plt.xlabel("x, m")
            plt.ylabel("Amplitude")
            plt.legend(loc="lower left")
            plt.grid()
            plt.show()
            return True
        except:
            print("An exception occurred ")
        else:
            return False

    def plot_number_dof_per_wavelength(self):
        try:
            #Number of points per wavelength
            plt.plot(self.x,self.nder,label="Numerical derivative, 2 points",marker="+",color="blue")
            plt.title(f"First derivative, error = {self.rms} $n_\lambda$ = {self.l/self.dx}")
            plt.xlabel("x, m")
            plt.ylabel("Amplitude")
            plt.legend(loc="lower left")
            plt.xlim((self.xmax/2-self.l,self.xmax/2+self.l))
            plt.grid()
            plt.show()
            return True
        except:
            print("An exception occurred ")
        else:
            return False
