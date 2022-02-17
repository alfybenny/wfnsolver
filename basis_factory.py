# -*- coding: utf-8 -*-
"""
Created on Sat Dec 25 16:26:21 2021

@author: alfys
"""
import basis

def create(basis_type, x, n):
    """
    Creates an object of a `prob_rule` class based on the given type.

    :param prob_rule_type: The type of the probability rule. Currently supported "arrhenius", "fret", "uniform", and "marcus".
    """
    
    if basis_type.lower() == 'harmonic':
        return basis.harm_pot(x, n)
    elif basis_type.lower() == '1d_box':
        return basis.particle_box(x, n)
    elif basis_type.lower() == 'num_harmonic':
        return basis.num_harm_pot(x, n)
    else:
        raise ValueError(format)
    