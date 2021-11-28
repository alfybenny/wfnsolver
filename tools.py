# -*- coding: utf-8 -*-
"""
Created on Fri Nov 26 21:51:38 2021

@author: alfys
"""
import numpy as np
import basis


# WARNING: make which basis to use choice of user
# def get_coeff(x, psi_x, phi_i, n, L):
def get_coeff(x, psi_x, n):
    cn = np.zeros(len(n), dtype = complex)
    dx = x[1]-x[0]
    
    for i in range (0, len(cn)):
        
        base = 0+0j
        # phi_n = basis.particle_box(x, n[i])
        phi_n = basis.harm_pot(x, n[i])
        phi_n_x = phi_n.wfn()

        for j in range (0, len(x)):
            base += phi_n_x[j]*psi_x[j]*dx

        cn[i] = base

    return cn

def get_initial_wpkt(x, n, cn, psi_x):
        
    initial_wpkt = np.zeros(len(x), dtype=complex)
    
    for i in range(0, len(n)):
        # phi_n = basis.particle_box(x, n[i])
        phi_n = basis.harm_pot(x, n[i])
        
        initial_wpkt += cn[i] * phi_n.wfn() # * prop(n[i], L, t[j])
        # c1 += cn[i]
        
    return initial_wpkt
        
        