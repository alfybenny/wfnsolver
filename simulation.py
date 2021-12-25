# -*- coding: utf-8 -*-
"""
Created on Fri Nov 26 20:42:21 2021

@author: alfys
"""
import numpy as np
import w_packet
import tools
from matplotlib import pyplot as plt

# Fundamental parameters for the simulation------
L = 4 # x axis from 0 to L
x_resolution = 2000 # x axis resolution
x = np.linspace(0, L, x_resolution)

# Total time for the simulation------------------
t = np.linspace(0, 10, 100)

# Number of basis sets used to define wavepacket-
n = np.linspace(1, 150, 150)

# Initial wavepacket-----------------------------
sigma = 0.5 # Gaussian wavepacket 
x0 = 1.7 # Initial position of wavepacket
k0 = 0.5 # Momentum of wavepacket

# Define wavepacket =
## 1. Gaussian: w_packet.gauss
## 2. Lorenzian (To be implimented)

wavepacket_type = w_packet.gauss(sigma, x, x0, k0)
psi_x = wavepacket_type.construct_wpkt() # Construct the wavepacket


# Basis_set--------------------------------------
## 1. Particle in a box: particle_box
## 2. Harmonic oscillator: harm_pot

basis_type = 'harmonic'
#------------------------------------------------


cn = tools.get_coeff(x, psi_x, n, basis_type) # Get coefficients

initial_wpkt = tools.get_initial_wpkt(x, n, cn, psi_x, basis_type) # Get initial wavepacket
                                                       # in terms of the basis
                                                       
# Plot and compare the original wavepacket and the fitted one
plt.plot(x, np.real(psi_x))
plt.plot(x, np.real(initial_wpkt))
plt.savefig('initial_wavepacket.png')
plt.close()


# Propogate--------------------------------------
tools.propogate(x, n, cn, t, basis_type)


