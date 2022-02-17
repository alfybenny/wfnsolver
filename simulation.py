# -*- coding: utf-8 -*-
"""
Created on Fri Nov 26 20:42:21 2021

@author: alfys
"""
import numpy as np
import w_packet
import tools
from matplotlib import pyplot as plt
import basis_factory

# Fundamental parameters for the simulation------
# L = 4 # x axis from 0 to L
L = 4 # x axis from 0 to L
x_resolution = 2000 # x axis resolution
x = np.linspace(0, L, x_resolution)

# Total time for the simulation------------------
t = np.linspace(0, 0.02, 1)

# Number of basis sets used to define wavepacket-
n = np.linspace(1, 400, 400)

# Initial wavepacket-----------------------------
sigma = 0.03 # Gaussian wavepacket 
# x0 = 1.7 # Initial position of wavepacket
x0 = 1.5 # Initial position of wavepacket
k0 = 0 # Momentum of wavepacket

# Define wavepacket =
## 1. Gaussian: w_packet.gauss
## 2. Lorenzian (To be implimented)

wavepacket_type = w_packet.gauss(sigma, x, x0, k0)
psi_x = wavepacket_type.construct_wpkt() # Construct the wavepacket


# Basis_set--------------------------------------
## 1. Particle in a box: particle_box
## 2. Harmonic oscillator: harm_pot

# Adding basis_set from numerical methods--------
## 1. Harmonic oscillator: num_harm_pot
####: WARNING: For numerical methods, you have to first store the eigenvalues 
####  and eigenvectors. This is implimented in the subclass.

# basis_type = 'harmonic'
# basis_type = '1d_box'
basis_type = 'num_harmonic'

# Numerical method eigenstoring:
## For now already done in the basis.py file



#------------------------------------------------


cn = tools.get_coeff(x, psi_x, n, basis_type) # Get coefficients

initial_wpkt = tools.get_initial_wpkt(x, n, cn, psi_x, basis_type) # Get initial wavepacket
                                                       # in terms of the basis
                                                       
# Plot and compare the original wavepacket and the fitted one
for i in range(0, 10):

    phi_n = basis_factory.create(basis_type, x, i+1)
    phi_n_x = phi_n.wfn()
#plt.plot(x, np.real(psi_x))
#plt.plot(x, np.real(initial_wpkt))
    plt.plot(x, phi_n_x)
plt.savefig('initial_wavepacket.png')
plt.close()


# Propogate--------------------------------------
#y_array = tools.propogate(x, n, cn, t, basis_type)



#tools.animate(x, y_array, basis_type, n, 'test.gif')

# tools.print_to_file(x, y_array, 'data.txt')


