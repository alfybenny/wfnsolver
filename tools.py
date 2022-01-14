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
    
    y_array = []
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
            
            
        y_array.append(abs(store_wave_t))
        store_wave_t = np.zeros(len(x), dtype=complex) # Reflush
        
    return y_array
 
        
def animate(x, y, basis_type, n, out_filename):
    for i in range(0, len(y)):
        basis_pot = basis_factory.create(basis_type, x, n)
        pot = basis_pot.potential()
        plt.plot(x, pot)
        
        plt.plot(x, y[i])
        plt.fill_between(x, y[i], step='pre', alpha = 0.4, color = 'orange')
        # plt.plot(x, abs(norm_wave_t))
        plt.ylim([0,1])
        plt.savefig(str(i)+'.png')
        plt.close()

    images = []

    for filename in range(0, len(y)):
        images.append(imageio.imread(str(filename)+'.png')) 
    imageio.mimsave(out_filename, images)
    
def print_to_file(x, y, filename):
    f = open(filename, 'w')
    for i in range(0, len(y)):
        f.write('step\n')
        j = y[i]
        for k in range(0, len(x)):
            f.write("  "+str(x[k])+" "+str(j[k])+"\n")
    f.close()
        
def animate_2d(x, y, x_array, y_array):
    
    def plot_3d(a, b, c, name):
        
        X, Y = np.meshgrid(a, b)
        
        fig = plt.figure()
        ax = plt.axes(projection='3d')
        ax.contour3D(X, Y, c, 50, cmap='jet')
        ax.set_zlim(0, 0.2)
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel('z')
        plt.savefig(str(name)+'.png')
        plt.close()
        
    for i in range(0, len(x_array)):
        z_matrix = []
        for n in x_array[i]:
            row = []
            for m in y_array[i]:
                row.append(n*m)
            z_matrix.append(row)
            
        plot_3d(x, y, z_matrix, i)
    images = []

    for filename in range(0, len(y)):
        images.append(imageio.imread(str(filename)+'.png')) 
        imageio.mimsave('out_filename.gif', images)    