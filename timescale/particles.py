"""
particles.py
author: Wolf Cukier
Methods to generate a set of initial condition state vectors for integration for
my Venus JP
"""

import spiceypy as spice
import numpy as np
from timescale.constants import *
import json

with open("data/planets.json") as f:
    planet_data = json.load(f)
spice.furnsh( 'data/meta.tm' )


def uniform_sphere(n):
    """
    Returns n points approximately uniformly spaced on a unit sphere using the
    Fibbonacci sphere algorithm.

    Units
    n:            dimentionless
    return:       dimentionless (n,3)

    Equation Source: https://stackoverflow.com/a/26127012
    """

    points = []
    phi =  np.pi* (3 - np.sqrt(5))

    ys = np.linspace(-1, 1, n)

    for i, y in enumerate(ys):
        r = np.sqrt(1 - y**2)
        theta = i*phi

        points.append((r*np.cos(theta), y, r*np.sin(theta)))

    return np.array(points)

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



def initial_state(n, v_inf, planets = ["2", "3", "5"], source_id="2"):
    """
    Returns the state vectors for the specified planets followed by n uniformly
    spaced particles about a sphere 300 Radii of the planet indicated by
    source_id away from that planet moving at v_inf relative to that planet.
    By default uses Venus, Earth, and Jupiter with Venus as the source.
    Alternative sets of planets can be specified by passing a list of SPICE IDs
    to planets.  Planets arragement is determined by picking a random time
    between Jan-1-1850 and Jan-1-2150.  All vectors are in ECLIP_J2000.
    Source_id must be in the array of planets.

    Units
    n:                dimentionless
    vel:              m/s
    planets:          dimentionless (List of Strings of the digets 1-8)
    source_id         dimentionless (String of a digit 1-8)
    return:           m, m/s (n + len(planets), 6)
    return[:,:3]:     m
    return[:,3:]:     m/s
    """
    if source_id not in planets:
        print("source_id not in planet array")
        raise

    time = random_time(1850, 2150)
    states = np.zeros((n+len(planets), 6))

    source = []

    for i, p in enumerate(planets):
        [pos, lt] = spice.spkezr(p, time, 'ECLIPJ2000', 'NONE', 'SUN')
        pos *= 1000 # convert from km to m
        if (p == source_id): source = pos
        states[i,:] = pos

    points = uniform_sphere(n)
    states[len(planets):,:] = np.hstack(
        (points * 100 * planet_data[source_id]["radius"], points * v_inf))

    states[len(planets):,:] += source


    return states