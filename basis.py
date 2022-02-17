# -*- coding: utf-8 -*-
"""
Created on Fri Nov 26 20:50:15 2021

@author: alfys
"""
# CLASS==========================================
import numpy as np
import numpy.polynomial.hermite as Herm
import math
from scipy.linalg import eigh_tridiagonal
from abc import abstractmethod

L = 4 # x axis from 0 to L
x_resolution = 2000 # x axis resolution
x = np.linspace(0, L, x_resolution)



def build_basis(x): # dx will be provided from driver code
    
    def V(x): # Define potential
        return 500000*((x-2)**4-(x-2)**2-0.0*(x-2)) # WARNING: the 0.5 should be automatically changed to what ever user inputs
    
    dx = 1/2000

    d = 1/dx**2 + V(x)[1:-1]
    e = -1/(2*dx**2)*np.ones(len(d)-1)
    
    w, v = eigh_tridiagonal(d, e)
    
    return w, v.T

a, b = build_basis(x)
    # num_basis_store = v.T
    # num_energy_store = w

class basis_set:
    
    num_basis_store = b # For storing basis functions for numerical methods so that
                         # dont have to calculate it again and again
    num_energy_store = a
    
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

## ANALYTICAL METHODS==========================================================
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
        pot = (self.x-5)**2
        return pot
        
    
class harm_pot(basis_set):
        
    def hermite(self, b, n):
        n = int(n)
        xi = np.sqrt(1*1/1)*(b)
        herm_coeffs = np.zeros(n+1)
        herm_coeffs[n] = 1
        return Herm.hermval(xi, herm_coeffs)
    
    def wfn(self):
        xi = np.sqrt(1*1/1)*(self.x - 0.5)
        prefactor = 1./math.sqrt(2.**self.n * math.factorial(self.n)) * (1*1/(np.pi*1))**(0.25)
        psi_n = prefactor * np.exp(- xi**2) * self.hermite((self.x - 2), self.n)

        

        return psi_n
    
    def energy(self, n):
        En = self.n + 1/2.0
        return En
    
    def  potential(self):
        pot = (self.x-2)**2
        return pot
    
class morse(basis_set):
    # Source: https://scipython.com/blog/the-morse-oscillator/
    def wfn(self):
        pass
    
    def potential(self):
        pass

## NUMERICAL METHODS===========================================================
class num_harm_pot(basis_set):
    
    
    
    def V(self, x): # Define potential
        return (x-0.5)**2 # WARNING: the 0.5 should be automatically changed to what ever user inputs
    
    def build_basis(self): # dx will be provided from driver code
        
        dx = 1/2000
    
        d = 1/dx**2 + self.V(self.x)[1:-1]
        e = -1/(2*dx**2)*np.ones(len(d)-1)
        
        w, v = eigh_tridiagonal(d, e)
        
        self.num_basis_store = v.T
        self.num_energy_store = w
        
        # return w, v.T
        
    
    def wfn(self):
        
        # dx = 1/2000
        
        # w, v = self.build_basis(dx, self.x)
        
        # M = v[int(self.n)]
        # print(self.n)
        M = self.num_basis_store[int(self.n-1)]
          # Have to add the first and last zeros
        M1 = np.append(M, 0.0) 
        M2 = np.insert(M1, 0, 0.0)
        M3 = M2*10
        
        return M3
    
    def potential(self):
        
        pot =((x-2)**4)-((x-2)**2)-(0.1*(x-2))
        
        return pot
    
    def energy(self, n):
                
        return self.num_energy_store[int(self.n-1)]
