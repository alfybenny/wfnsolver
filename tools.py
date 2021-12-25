# -*- coding: utf-8 -*-
"""
Created on Fri Nov 26 21:51:38 2021

@author: alfys
"""
import numpy as np
import basis
import basis_factory
from matplotlib import pyplot as plt
import imageio


# WARNING: make which basis to use choice of user
# def get_coeff(x, psi_x, phi_i, n, L):
    
def get_coeff(x, psi_x, n, basis_type):
    cn = np.zeros(len(n), dtype = complex)
    dx = x[1]-x[0]
    
    for i in range (0, len(cn)):
        
        base = 0+0j
        
        phi_n = basis_factory.create(basis_type, x, n[i])
        #######################################################################
        # phi_n = basis.particle_box(x, n[i])
        # phi_n = basis.harm_pot(x, n[i])
        #######################################################################
        phi_n_x = phi_n.wfn()

        for j in range (0, len(x)):
            base += phi_n_x[j]*psi_x[j]*dx

        cn[i] = base

    return cn

def get_initial_wpkt(x, n, cn, psi_x, basis_type):
        
    initial_wpkt = np.zeros(len(x), dtype=complex)
    
    for i in range(0, len(n)):
        phi_n = basis_factory.create(basis_type, x, n[i])
        #######################################################################
        # phi_n = basis.particle_box(x, n[i])
        # phi_n = basis.harm_pot(x, n[i])
        #######################################################################
        
        initial_wpkt += cn[i] * phi_n.wfn() # * prop(n[i], L, t[j])
        # c1 += cn[i]
        
    return initial_wpkt

def propogate(x, n, cn, t, basis_type):
    
    def prop(x, n, t):
        basis_energy = basis_factory.create(basis_type, x, n)
        #######################################################################
        # basis_type = basis.particle_box(x, n)
        # basis_type = basis.harm_pot(x, n)
        #######################################################################
        E = basis_energy.energy(n)
        ci = 0.0 +1j
        psi_n_t = np.exp(-1*ci*E*t)
        return psi_n_t
        
    store_wave_t = np.zeros(len(x), dtype=complex)
    
    for j in range(0, len(t)):
        for i in range(0, len(cn)):
            phi_n = basis_factory.create(basis_type, x, n[i])
            #######################################################################
            # phi_n = basis.particle_box(x, n[i])
            # phi_n = basis.harm_pot(x, n[i])
            #######################################################################
            store_wave_t += cn[i] * phi_n.wfn() * prop(x, n[i], t[j])
            
        norm_constant = np.linalg.norm(store_wave_t)
        # norm_wave_t = store_wave_t / norm_constant
        
        # Plotting the potential energy
        basis_pot = basis_factory.create(basis_type, x, n)
        pot = basis_pot.potential()
        plt.plot(x, pot)
        
        plt.plot(x, abs(store_wave_t))
        # plt.plot(x, abs(norm_wave_t))
        plt.ylim([0,1])
        plt.savefig(str(j)+'.png')
        plt.close()
        store_wave_t = np.zeros(len(x), dtype=complex) # Reflush
    
    images = []

    for filename in range(0, len(t)):
        images.append(imageio.imread(str(filename)+'.png')) 
    imageio.mimsave('wave1.gif', images)
        # phi_n = basis.particle_box(x, n[i])
        
        
        # initial_wpkt += cn[i] * phi_n.wfn() # * prop(n[i], L, t[j])
        # c1 += cn[i]
        
    # return initial_wpkt
        

    
    