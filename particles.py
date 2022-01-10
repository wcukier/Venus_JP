"""
particles.py
author: Wolf Cukier
Methods to generate a set of initial condition state vectors for integration for
my Venus JP
"""

import spiceypy as spice
import numpy as np
from constants import *


def random_direction():
    """
    Returns a 3-D unit vector pointing in a random direction
    
    Units
    return:      dimentionless
    """
    
    theta = np.random.random() * np.pi * 2
    phi = np.random.random() * np.pi
    
    x = np.cos(theta)*np.sin(phi)
    y = np.sin(theta)*np.sin(phi)
    z = np.cos(phi)
    
    return np.array([x, y, z])


def random_time(start, end):
    """
    Takes a start year, start, and an end year, end, and returns a random 
    ephemeris time between them.  B.C should be a negative number.
    
    Units
    start:        calendar year
    end:          calendar year
    return:       ephemeris time
    """
    
    if start < 0: start = str(start) + " B.C."
    else: start = str(start) + " A.D."
    
    if end < 0: end = str(end) + " B.C."
    else: end = str(end) + " A.D." 

    
    start = spice.str2et(f"{start} Jan 1")
    end = spice.str2et(f"{end} Jan 1")
    
    return start + int(np.random.random() * (end-start))
    
    
    
    