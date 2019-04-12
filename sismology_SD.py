import numpy as np
from math import *
import matplotlib.pyplot as plt

class Sismology_SD():
    """
    Example for second derivative using finite difference method (by simology course)

    https://www.coursera.org/learn/computers-waves-simulations/

    """
    def __init__(self,xmax,nx):
        #initial variables
        self.xmax = xmax                          #physical domain
        self.nx = nx                              #number of samples
        self.dx = self.grid_dx()                  #grid increment
        self.x = np.linspace(0,self.xmax,self.nx) #space coordinates
        self.rms = 0.

        # Initiation of numerical and analytical derivatives
        self.nder3 = np.zeros(self.nx)             #numerical derivative
        self.nder5 = np.zeros(self.nx)             #numerical derivative
        self.ader = np.zeros(self.nx)             #analytical derivative

        #Initiation of Gaussian function
        self.x0 = self.xmax/2                     # Center of Gaussian function x0 (m)
        self.a=.25                                # exponent of Gaussian function
        self.f = self.analytical_solution()       # Gaussian function

    def grid_dx(self):
        """
        this function create a uniform distribution
        """
        return self.xmax/(self.nx-1)

    def analytical_solution(self):
        """
        Initialization of Gaussian function
        """
        return (1./sqrt(2*pi*self.a))*np.exp(-(((self.x-self.x0)**2)/(2*self.a)))

    def analytical_second_derivative_solution(self):
        """
        analytical second derivative of given function
        """
        self.ader = 1./sqrt(2*pi*self.a)*((self.x-self.x0)**2/self.a**2 -1/self.a)*np.exp(-1/(2*self.a)*(self.x-self.x0)**2)
        #exclude boundaries
        self.ader[0] = 0.0
        self.ader[self.nx-2] = 0.0


    def second_derivative_5pts(self):
        """
        second derivative with 5 points
        """

        #numerical second derivative of a given Function
        for i in range (1, self.nx-2):
            self.nder5[i] = (-1./12 * self.f[i - 2] + 4./3  * self.f[i - 1] - 5./2 * self.f[i] \
                       +4./3  * self.f[i + 1] - 1./12  * self.f[i + 2]) / self.dx ** 2

    def second_derivative_3pts(self):
        """
        second derivative with 3 points
        """

        #numerical second derivative of a given Function
        for i in range(1,self.nx-1):
            self.nder3[i]=(self.f[i+1] - 2*self.f[i] + self.f[i-1])/(self.dx**2)

    def error_5pts(self):
        """
        class error calculation: max norm
        """
        self.rms = self.rms *0
        self.rms = np.sqrt(np.mean(self.nder5-self.ader)**2)

    def error_3pts(self):
        """
        class error calculation: max norm
        """
        self.rms = self.rms *0
        self.rms = np.sqrt(np.mean(self.nder3-self.ader)**2)


    def plot_analytic_wave(self):
        try:
            plt.plot(self.x,self.f)
            plt.title("Gaussian function")
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
            plt.figure(figsize=(10,6))
            plt.plot(self.x,self.nder5,label="Numerical second derivative, 5 points",marker="+",color="violet")
            plt.plot(self.x,self.nder3,label="Numerical second derivative, 3 points",marker="o")
            plt.plot(self.x,self.ader,label="Analytical derivative",lw=2,ls="--")
            plt.plot(self.x,self.nder5-self.ader,label="Diference",lw=2,ls=":")
            plt.title(f"Second derivative, Err (rms) = {self.rms}")
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
