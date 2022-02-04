"""
collisions.py
author: Wolf Cukier
Methods to caclculate the collision probability of a particle with a target 
body with arbitrary eccentricity
"""

import numpy as np
from timescale.constants import *
from scipy.integrate import quad


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
    # if (A < 1): return 0
   
    result = np.sqrt((A**2 * e**2 - (A - 1)**2) / (A**2 * (1 - e**2)))
    # if (np.isnan(result)): result = 0
    
    return result

def rel_vel(a, e, i, rho, theta, target):
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

    a_0 = target["semi_major"]
    e_0 = target["eccentricity"]

    A = a/rho
    A_0 = a_0/rho
    cot_alph = cot_alpha(A,e)
    cot_alph_0 = cot_alpha(A_0, e_0)
    if (theta < np.pi): sgn = 1
    else: sgn = -1


    Usq = G * MASS_SUN/(rho**2) * (a * (2*A - 1)/(A**2)
                                   + a_0 * (2*A_0 - 1)/(A_0**2)
                                   - 2* np.sqrt(a*a_0*(1-e**2)*(1-e_0**2))
                                   * (np.cos(i) + sgn*cot_alph*cot_alph_0))
    return np.sqrt(np.abs(Usq))


def radius(theta, target):
    """
    Returns the radial distance from the Sun of the target for a particular 
    value of theta.
    
    Inputs
    target:     dict
    theta:      radians
    
    Output
    rho:        m

    Source: Wetherill (1969, eqn. 25)
    """
    
    a_0 = target["semi_major"]
    e_0 = target["eccentricity"]
    
    return a_0 * (1 - e_0**2) / (1 + e_0 * np.cos(theta))
    



def prob_per_theta(theta, a, e, i, target):
    """
    Takes the semi-major axis, a, eccentricity, e, and orbital inclination, i,
    of a particle and returns the probability per unit time per unit angle that
    the particle will collide with the target specified in config.  This
    measurement is only valid over long time scales.
    
    Units
    a:         m
    e:         dimentionless
    i:         radians
    target:    dictionary
    return:    s^-1

    Equation Source: Wetherill (1969)
    """

    tau = target["radius"]
    a_0 = target["semi_major"]
    e_0 = target["eccentricity"]
    rho = radius(theta, target)
    A = a/rho
    abs_cot_theta = np.abs(cot_alpha(A,e))

    A_0 = a_0/rho
    U = rel_vel(a, e, i, rho, theta, target)


    P_F =  (((tau**2) * U * rho)
            / (8 * (np.pi**2) * np.sin(i) * abs_cot_theta * (a_0**2)
               * np.sqrt(1 - e_0**2) * (a**2) * np.sqrt(1 - e**2)))
    if np.isnan(P_F): return 0
    if P_F < 0: return 0
    if P_F >= 1: return 0

    return P_F




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
    i:         radians
    target:    dictionary
    return:    s^-1
    
    Equation Source: Wetherill (1969)
    """
    # return prob_per_theta(a, a, e, i, target)
    p, _ = quad(prob_per_theta, 0, 2*np.pi, args=(a, e, i, target))
    return p
    
    


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
