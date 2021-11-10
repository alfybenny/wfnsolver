import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation

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

x0 = 200.0 # Initial position of wavepacket
k0 = 0.4 # Average Momentum
sig = 15.0 # Spread of initial wavepacket

n = np.linspace(1, 100, 100) # Number of quantized basis functions used for the simulation

# Generation of wavepacket=================================
wavepacket = gaus_fun(x, x0, sig, k0)

# Calculation of fourier coefficients----------------------
cn = FA(x, wavepacket, n, L)

# base to add each coefficients one by one
psi_total = np.zeros(len(x), dtype=complex)

# Summation of n basis functions with fourier coefficients
for i in range(0, len(cn)):
    psi_total += cn[i] * phi(x, n[i], L)

# Time evolution of wavepacket=============================
psi_total_time = np.zeros(len(x), dtype=complex) # Empty array
sig_cn = np.zeros(len(x), dtype=complex)
t = np.linspace(0, 10000, 500) # Total time and timesteps

for j in range(0, len(t)):
    for i in range(0, len(cn)):
        psi_total_time += cn[i] * phi(x, n[i], L) * prop(n[i], L, t[j])
        sig_cn += cn[i]*np.conj(cn[i])
        
    norm_psi_total_time = psi_total_time/sig_cn
    # Output-----------------------------------------------
    plt.ylim(-1,1)
    plt.plot(x, np.real(norm_psi_total_time), 'red')
    plt.savefig('dyna_'+str(j)+'.png')
    plt.close()
    psi_total_time = np.zeros(len(x), dtype=complex) # Restarting dynamic array for t + dt
    sig_cn = np.zeros(len(x), dtype=complex)

# Plot fourier coefficients--------------------------------
plt.close()
plt.plot(n, np.real(cn))
plt.savefig('coeff.png')

#END=======================================================

