# -*- coding: utf-8 -*-
"""
Created on Fri Nov 26 20:50:15 2021

@author: alfys
"""
# CLASS==========================================
import numpy as np
import numpy.polynomial.hermite as Herm
import math
from abc import abstractmethod

class basis_set:
    
    def __init__(self, x, n):
        self.x = x
        # self.L = L
        self.n = n
    
    @abstractmethod
    def wfn():
        pass

    @abstractmethod
    def energy(self, m):
        pass
    
    @abstractmethod
    def potential():
        pass

# SUBCLASS=======================================
    
class particle_box(basis_set):
        
    def L(self, x):
        return x[len(x)-1]
    # L = self.x[len(self.x)-1]
    def wfn(self):
        # self.L = self.x[len(self.x)-1]
        phi_n = np.sqrt(2.0/self.L(self.x))*np.sin((self.n*np.pi*self.x)/self.L(self.x))
        return phi_n
    
    def energy(self, n):
        En = (self.n**2 * np.pi**2)/(2*self.L(self.x)**2)
        return En
    
    def potential(self):
        pass
        
    
class harm_pot(basis_set):
        
    def hermite(self, x, n):
        n = int(n)
        xi = np.sqrt(1*1/1)*x
        herm_coeffs = np.zeros(n+1)
        herm_coeffs[n] = 1
        return Herm.hermval(xi, herm_coeffs)
    
    def wfn(self):
        xi = np.sqrt(1*1/1)*(self.x - 5)
        prefactor = 1./math.sqrt(2.**self.n * math.factorial(self.n)) * (1*1/(np.pi*1))**(0.25)
        psi_n = prefactor * np.exp(- xi**2) * self.hermite((self.x - 5), self.n)
        # norm_psi_n = np.linalg.norm(psi_n)
        # final_psi_n = psi_n / norm_psi_n
        # return final_psi_n
        return psi_n
    
    def energy(self, n):
        En = self.n + 1/2.0
        return En
    
    def  potential(self):
        pot = (self.x-5)**2
        return pot