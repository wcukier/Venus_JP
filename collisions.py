"""
collisions.py
author: Wolf Cukier
Methods to caclculate the collision probability of a particle with earth
"""

import numpy as np
from .constants import *

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
    return np.sqrt((A**2 * e**2 - (A - 1)**2) / (A**2 * (1 - e**2)))

def rel_vel(A, e, i):
    """
    Takes the scaled semimajor axis, A, eccentricity, e, and orbital 
    inclination, i, of a particle and returns the relative velocity between the 
    particle and Earth at the point of collision. 
    
    Units
    A:        dimentionless
    e:        dimentionless
    i:        degrees
    return:   m/s
    
    Equation Source: Öpik (1951)
    """
    return (np.sqrt(G * MASS_SUN / SEMI_MAJOR_EARTH) * 
            (3 + 1/A - 2 * np.sqrt(A * (1-e) **2) * np.cos(i)))


def prob_per_time(a, e, i):
    """
    Takes the semi-major axis, a, eccentricity, e, and orbital inclination, i, 
    of a particle and returns the probability per unit time that the particle
    will collide with earth.  This measurement is only valid over long time 
    scales.  This assumes that the radius of Earth is much larger than the 
    particle radius.
    
    Units
    a:         m
    e:         dimentionless
    i:         degrees
    return:    s^-1
    
    Equation Source: Wetherill (1969)
    """
    
    A = a/SEMI_MAJOR_EARTH
    U = rel_vel(A, e, i)
    abs_cot_alpha = np.abs(cot_alpha(A, e))
    
    return (RADIUS_EARTH**2 * U) / (2 * np.pi * np.sin(i) * SEMI_MAJOR_EARTH
                                    * a**2 * np.sqrt(1-e**2) * abs_cot_alpha)
    
    
def collision_probability(a, e, i, t):
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
    
    P_F = prob_per_time(a, e, i)
    
    return 1 - np.exp(-P_F * t)
    