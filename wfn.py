# Author: Alfy Benny

import numpy as np
import math
from sympy import *
import argparse
import os.path

V, x, L, n, i, h, m = symbols('V x L n i h m')

L = 1
n = 4
#V = x**2

def psi(i, x):
    return sqrt(2/L)*sin(i*pi*(x+L/2)/L)

def H_phi(i,x):
    return -(h/2*m)*diff(psi(i,x),x,x)+(V*psi(i,x))

def dotproduct(i, j):
    return integrate(psi(i, x)*H_phi(j, x), (x, -1, 1))

range = [1, 2, 3];


def f(i, j):
    return dotproduct(i+1, j+1)

def basis(i, j):
    return psi(i+1, j+1)


M = Matrix(n, n, f)

B = Matrix(n, n, basis)

print(M)
#print(dotproduct(1,1))
E = M.eigenvals()
V = M.eigenvects()

wfn = V*B


print(wfn)

