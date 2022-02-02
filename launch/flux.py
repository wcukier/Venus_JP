"""
flux.py
author: Wolf Cukier

A set of equations that return the asteroid flux at a paticular time above a
certain diameter
"""

import numpy as np

def robbins_1km(t, ignore_unconstrained=False):
    """
    Returns the asteroid flux for D > 1 km at a time t years ago using the
    model in Robbins (2013).  If ignore_constrained is false or not set, will 
    throw an error if t any t is older than 3.92 Ga
    
    Inputs
    t:                          Gyr
    ignore_unconstrained:       Boolean
    
    Output
    flux                        km^-2 yr^-1

    Source: Robbins (2013) (Take the derivitive as per Neukam 2001)
    """
    
    
    alpha = 7.26e-41
    beta = 22.6
    gamma = 9.49e-4
    delta = 1.88e-4
    
    if (np.any(t > 3.98) and not ignore_unconstrained): raise
    
    return alpha * beta * np.exp(beta * t) + gamma + 2 * delta * t

def nesvorny_10km(t, model = 0):
    """
    Returns the asteroid flux for D > 10 km incident on earth at a time t years
    ago using the models in Nesvorny (2017).  Model 0 is the REF model, 1 is the
    CASE1B model.

    Inputs:
    t:                          Gyr
    model:                      0 or 1

    Output
    flux                       yr^-1 (at earth)

    Source: Nesvorny (2017)
    """
    if (model == 0):
        tau = [.037, .16, 1.5]
        F = [0.5e3, .31e3, .0027e3]
    elif (model == 1):
        tau = [.027, .16, 100]
        F = [2.9e3, .6e3, .0006e3]
    else:
        raise
    flux = 0
    for i in range(3):
        flux += F[i] * np.exp(-(4.5 - t)/tau[i])

    return flux


