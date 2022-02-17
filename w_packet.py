
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 26 20:26:35 2021

@author: alfys
"""
import numpy as np
from abc import abstractmethod

# All the wavepackets implimented in this simulation will be under the following class
class wave_packet:
    
    # attributes
    def __init__(self, sigma, x, x0, k0):
        self.sigma = sigma
        self.x = x
        self.x0 = x0
        self.k0 = k0
    
    # virtual methods
    @abstractmethod
    def construct_wpkt(self, sigma, x, x0, k0):
        pass
  
# SUBCLASS=======================================
## 1. Gaussian wavepacket------------------------
class gauss(wave_packet):
    
    def construct_wpkt(self):
        ci = 0. + 1j
        pre = 1.0/(self.sigma*np.sqrt(2.0*np.pi))
        psi_x = pre*np.exp(-0.5*((self.x-self.x0)/self.sigma)**2)*np.exp(ci*self.k0*self.x)
        
        psi_x_N = psi_x/np.linalg.norm(psi_x)
        
        return psi_x

## 2. lorentian wavepacket-----------------------
class lorentz_wpkt(wave_packet):
    pass
    