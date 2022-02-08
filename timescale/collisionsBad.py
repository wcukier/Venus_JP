"""
collisions.py
author: Wolf Cukier
Methods to caclculate the collision probability of a particle with earth
"""

import numpy as np
from timescale.constants import *
import json


def cot_alpha(A, e):
    """
    Takes the scaled semimajor axis, A, and eccentricity, e, of the particle
    and returns the angle between the trajectory of the particle and the
    radius vector from the Sun.
    
    Units and Type
    A:         dimentionless
    e:         dimentionless
    return:    dimentionless
    
    Equation Source: Wetherill (1969)
    """
    # if (A < 1): return 0 #This is wrong
   
    result = np.sqrt((A**2 * e**2 - (A - 1)**2) / (A**2 * (1 - e**2)))
    # if (np.isnan(result)): result = 0 #This is wrong
    
    return result

def rel_vel(A, e, i, target):
    """
    Takes the scaled semimajor axis, A, eccentricity, e, and orbital
    inclination, i, of a particle and returns the relative velocity between the
    particle and the target at the point of collision.
    
    Units
    A:        dimentionless
    e:        dimentionless
    i:        degrees
    return:   m/s
    
    Equation Source: Öpik (1951)
    """

    # if (A < 1): return np.nan #This is wrong
    return (np.sqrt(G * MASS_SUN / target["semi_major"]) *
            (3 - 1/A - 2 * np.sqrt(A * (1-e**2)) * np.cos(i))) #sign error and (1-e**2)





def prob_per_time(a, e, i, target):
    """
    Takes the semi-major axis, a, eccentricity, e, and orbital inclination, i, 
    of a particle and returns the probability per unit time that the particle
    will collide with the target specified in config.  This measurement is only
    valid over long time scales.  This assumes that the radius of the target is
    much larger than the particle radius and that the target has eccentricity
    and orbital inclination of 0.
    
    Units
    a:         m
    e:         dimentionless
    i:         degrees
    target:    dictionary
    return:    s^-1
    
    Equation Source: Wetherill (1969)
    """
    
    A = a/target["semi_major"]
    U = rel_vel(A, e, i, target)
    abs_cot_alpha = np.abs(cot_alpha(A, e))
    
    # if (A < 1): return 0 #This is wrong
    return (target["radius"]**2 * U) / (2 * np.pi**2 * np.sin(i) # missing a pi
                                        * target["semi_major"] * a**2 *
                                        np.sqrt(1-e**2) * abs_cot_alpha)
    
    


def collision_probability(a, e, i, t, target):
    """
    Takes the semi-major axis, a, eccentricity, e, and orbital inclination, i,
    of a particle and returns the probability of that particle hitting Earth
    over some long time scale, t.
    
    Units
    a:         m
    e:         dimentionless
    i:         degrees
    t:         s
    return:    dimentionless
    """
    P_F = prob_per_time(a, e, i, target)
    
    return 1 - np.exp(-P_F * t)
