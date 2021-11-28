import numpy as np
from mpl_toolkits import mplot3d
from matplotlib import pyplot as plt
import imageio

# FUNCTIONS================================================
# Gaussian wavepacket generator----------------------------
def gaus_fun(x, x0, sig, k0):

    ci = 0. + 1j
    pre = 1.0/(sig*np.sqrt(2.0*np.pi))
    psi_x = pre*np.exp(-0.5*((x-x0)/sig)**2)*np.exp(ci*k0*x)
    return psi_x

# Defining a function of complete basis set----------------
def phi(x, n, L):
    phi_n = np.sqrt(2.0/L)*np.sin((n*np.pi*x)/L)
    return phi_n

# Fourier transform function-------------------------------
def FA(x, psi_x, n, L):
    cn = np.zeros(len(n), dtype = complex)
    dx = x[1]-x[0]
    for i in range (0, len(cn)):

        base = 0+0j

        phi_i = phi(x, n[i], L)

        for j in range (0, len(x)):
            base = base + phi_i[j]*psi_x[j]*dx

        cn[i] = base

    return cn

# Energy for nth state-------------------------------------
def psi_En(n, L):
    En = (n**2 * np.pi**2)/(2*L**2)

    return En

# Time propogator------------------------------------------
def prop(n, L, t):
    E = psi_En(n, L)
    ci = 0.+1j
    psi_n_t = np.exp(-1*ci*E*t)
    ### Write code here to define phi_n_t
    return psi_n_t

# Normalize------------------------------------------------
def Normalize(c):
    sm=0.
    for i in range(0, len(c)):
      sm=sm+np.conj(c[i])*c[i]

    cn = c/sm
    return cn


# Initial values===========================================
L = 500 # Length of the 1D box
x = np.linspace(0, L, 2000) # Resolution of x axis
y = np.linspace(0, L, 2000) # Resolution of y axis
x0 = 100.0 # Initial position of wavepacket
y0 = 200.0 
k0 = -0.4 # Average Momentum
k0y = -0.5
sig = 18.0 # Spread of initial wavepacket

n = np.linspace(1, 100, 100) # Number of quantized basis functions used for the simulation

# Generation of wavepacket=================================
wavepacket = gaus_fun(x, x0, sig, k0)
wavepacket_y = gaus_fun(y, y0, sig, k0y)

# Calculation of fourier coefficients----------------------
cn = FA(x, wavepacket, n, L)
cn_y = FA(y, wavepacket_y, n, L)

# base to add each coefficients one by one
psi_total = np.zeros(len(x), dtype=complex)
psi_total_y = np.zeros(len(y), dtype=complex)

# Summation of n basis functions with fourier coefficients
for i in range(0, len(cn)):
    psi_total += cn[i] * phi(x, n[i], L)

for i in range(0, len(cn_y)):
    psi_total_y += cn_y[i] * phi(y, n[i], L)

# Time evolution of wavepacket=============================
psi_total_time = np.zeros(len(x), dtype=complex) # Empty array
psi_total_time_y = np.zeros(len(y), dtype=complex)

psi_total_time_xy = np.zeros((len(x),len(y)), dtype=complex)

#sig_cn = np.zeros(len(x), dtype=complex)
#sig_cn_y = np.zeros(len(y), dtype=complex)

c1 = 0
c2 = 0
t = np.linspace(0, 10000, 100) # Total time and timesteps

for j in range(0, len(t)):
    for i in range(0, len(cn)):
        psi_total_time += cn[i] * phi(x, n[i], L) * prop(n[i], L, t[j])
        c1 += cn[i]
    for k in range(0, len(cn_y)):
        psi_total_time_y += cn_y[k] * phi(y, n[k], L) * prop(n[k], L, t[j])
        c2 += cn_y[k]
        

    for a in range(0, len(x)):
        for b in range(0, len(y)):
            psi_total_time_xy[a][b] = (psi_total_time[a]/abs(c1)) * (psi_total_time_y[b]/abs(c2))
        
    
    norm_psi_total_time_xy = (abs(psi_total_time_xy))**2
    
    # Output-----------------------------------------------
    
    #plt.matshow(norm_psi_total_time_xy)
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    
    ax.plot_surface(x, y, norm_psi_total_time_xy, rstride=1, cstride=1, cmap='viridis', edgecolor='none')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z');


    plt.savefig(str(j)+'.png')
    plt.close()
    psi_total_time_xy = np.zeros((len(x),len(y)), dtype=complex)
    psi_total_time = np.zeros(len(x), dtype=complex)
    psi_total_time_y = np.zeros(len(y), dtype=complex)

images = []

for filename in range(0, len(t)):
    images.append(imageio.imread(str(filename)+'.png'))
imageio.mimsave('wave.gif', images)


# Plot fourier coefficients--------------------------------
plt.close()
plt.plot(n, np.real(cn))
plt.savefig('coeff.png')

#END=======================================================

