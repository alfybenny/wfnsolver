# -*- coding: utf-8 -*-
"""
Created on Thu Jan 13 19:29:59 2022

@author: alfys
"""

import numpy as np
import w_packet
import tools
from matplotlib import pyplot as plt

# Fundamental parameters for the simulation------
L = [10, 10]

resolution = [2000, 2000]

x = np.linspace(0, L[0], resolution[0])
y = np.linspace(0, L[1], resolution[1])

# Total time for the simulation------------------
t = np.linspace(0, 10, 50)

# Number of basis sets used to define wavepacket-
n = np.linspace(1, 150, 150)

# Initial wavepacket-----------------------------
sigma = [1, 1]
ini_pos = [5, 5]
k0 = [5, 5]

# Define wavepacket =
## 1. Gaussian: w_packet.gauss
## 2. Lorenzian (To be implimented)

wavepacketx_type = w_packet.gauss(sigma[0], x, ini_pos[0], k0[0])
wavepackety_type = w_packet.gauss(sigma[1], x, ini_pos[1], k0[1])
psi_x = wavepacketx_type.construct_wpkt() # Construct the wavepacket
psi_y = wavepackety_type.construct_wpkt()

# Basis_set--------------------------------------
## 1. Particle in a box: particle_box
## 2. Harmonic oscillator: harm_pot

# basis_type = 'harmonic'
basis_type = '1d_box'
#------------------------------------------------


cn_x = tools.get_coeff(x, psi_x, n, basis_type) # Get coefficients

initial_wpkt_x = tools.get_initial_wpkt(x, n, cn_x, psi_x, basis_type) # Get initial wavepacket
                                                       # in terms of the basis

cn_y = tools.get_coeff(y, psi_y, n, basis_type) # Get coefficients

initial_wpkt_y = tools.get_initial_wpkt(y, n, cn_y, psi_y, basis_type)                                                       
                                                       
# Plot and compare the original wavepacket and the fitted one
'''
plt.plot(x, np.real(psi_x))
plt.plot(x, np.real(initial_wpkt))
plt.savefig('initial_wavepacket.png')
plt.close()
'''

# Propogate--------------------------------------
x_array = tools.propogate(x, n, cn_x, t, basis_type)
y_array = tools.propogate(x, n, cn_y, t, basis_type)



# tools.animate(x, y_array, basis_type, n, 'test.gif')

tools.animate_2d(x, y, x_array, y_array)


