# -*- coding: utf-8 -*-
"""
Created on Fri Nov 26 20:42:21 2021

@author: alfys
"""
import numpy as np
import w_packet
import basis
import tools
from matplotlib import pyplot as plt

# Initial wavepacket-----------------------------
L = 10
x = np.linspace(0, L, 2000)
sigma = 0.3
x0 = 4.0
k0 = -10

wavepacket_type = w_packet.gauss(sigma, x, x0, k0)
psi_x = wavepacket_type.construct_wpkt()
# print(psi_x.construct_wpkt())

# Basis_set--------------------------------------

n = np.linspace(1, 100, 100)

basis_type = 'particle_box'

cn = tools.get_coeff(x, psi_x, n)

initial_wpkt = tools.get_initial_wpkt(x, n, cn, psi_x)

plt.plot(x, np.real(psi_x))
plt.plot(x, np.real(initial_wpkt))
plt.show()
plt.close()

# test = basis.harm_pot(x, 1000)
# test1 = test.wfn()

# plt.plot(x, test1)
# plt.show()
