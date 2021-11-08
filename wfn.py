# Author: Alfy Benny

import numpy as np
import math
from sympy import *
import argparse
import os.path

V, x, L, n, i, h, m = symbols('V x L n i h m')

L = 1

def psi(i, x):
    return sqrt(2/L)*sin(i*pi*(x+L/2)/L)

#psi = sqrt(2/L)*sin(i*pi*(x+L/2)/L)
def H_phi(i,x):
    return -(h/2*m)*diff(psi(i,x),x,x)-(V*psi(i,x))

#Hphi = -(h/2*m)*diff(psi,x,x) - V*psi



print(H_phi(2,x))


