"""
collisions.py
author: Wolf Cukier
Methods to caclculate the collision probability of a particle with earth
"""

import numpy as np
from timescale.constants import G, MASS_SUN
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
    # if ((A**2 * e**2) < ((A-1)**2)): return np.nan
    # result = np.sqrt((A**2 * e**2 - (A - 1)**2) / (A**2 * (1 - e**2)))

    rsq = ((A**2) * (e**2) - (A - 1)**2 )/ ((A**2) * (1-e**2))
    # print(f"cot: {result}")
    return np.sqrt(rsq)

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

    # if (A < 1): return np.nan


    # result = (np.sqrt((G * MASS_SUN / target["semi_major"]) *
    #         (3 - 1/A - 2 * np.sqrt(A * (1 - (e**2))) * np.cos(i))))

    # result = (np.sqrt((G * MASS_SUN / target["semi_major"]) *
    #         (3 - 1/A - 2 * np.sqrt(A * (1 - (e**2))) * np.cos(i))))

    Usq = (G * MASS_SUN / target["semi_major"]) * (3 - 1/A - 2 *np.cos(i) *np.sqrt(A * (1 - e**2)))

    # print(f'U:{result}')
    return np.sqrt(Usq)

# def rel_vel(A, e, i, target):
#     """
#     Takes the scaled semimajor axis, A, eccentricity, e, and orbital
#     inclination, i, of a particle and returns the relative velocity between the
#     particle and the target at the point of collision.

#     Units
#     A:        dimentionless
#     e:        dimentionless
#     i:        degrees
#     return:   m/s

#     Equation Source: Öpik (1951)
#     """

#     a_0 = target["semi_major"]
#     e_0 = 0
#     rho = a_0

#     # print(A)
#     # A_0 = a_0/a_0
#     # cot_alph = cot_alpha(A,e)
#     # cot_alph_0 = cot_alpha(A_0, e_0)
#     # sgn=1




#     Usq = np.abs(G * MASS_SUN/(rho) * ((2*A - 1)/(A)
#                                    + (1)
#                                    - (2* np.sqrt(A*(1-e**2)) * (np.cos(i)))))

#     print(f"U^2: {Usq}")
#     if (Usq < 0): return np.nan
#     return np.sqrt(Usq)




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
    if(np.isnan(a)):return 0
    A = a/target["semi_major"]
    # q = A*(1 - e)
    # Q = A*(1 + e)
    # if((q > 1) or (Q<1)): return 0
    # print(f"a: {a}, e: {e}, i: {i}")
    U = rel_vel(A, e, i, target)
    abs_cot_alpha = np.abs(cot_alpha(A, e))


    P_F = ((target["radius"]**2) * U) / (2 * np.pi**2 * np.sin(i) * target["semi_major"] * (a**2) * np.sqrt((1-(e**2))) * abs_cot_alpha)















    # if ((A*(1+e)) < 1): return 0
    # if ((A*(1-e)) > 1): return 0
    # print(type(U))
    # if (A < 1): return 0
    # P_F = (target["radius"]**2 * U) / (2 * np.pi**2 * np.sin(i)
    #                                     * target["semi_major"]
    #                                     * a**2
    #                                     * np.sqrt(1-e**2) * abs_cot_alpha)

    # if np.isnan(P_F): return 0
    # if P_F < 0: return 0
    # if P_F >= 1: return 0
    return P_F



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
    i = i + 1.57 * np.pi/180
    P_F = prob_per_time(a, e, i, target)
    return 1 - np.exp(-P_F * t)
