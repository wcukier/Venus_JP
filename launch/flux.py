"""
flux.py
author: Wolf Cukier

A set of equations that describe impactor flux on Venus
"""

from matplotlib.pyplot import sca
import numpy as np
import fluids
import json
with open("../data/planets.json") as f:
    planets = json.load(f)
import numba

VENUS_RADIUS = planets["2"]["radius"]
VENUS_MASS = planets["2"]["mass"]
rho_0 = fluids.ATMOSPHERE_1976(Z=0).rho


STEPHAN_BOLTZMAN = 5.67e-8 #W/(m^2 K^4)
GRAV_CONSTANT = 6.674e-11 # N m^2/kg^2

def N_m(T, m, alpha=4.09e-31, beta=22.6, gamma=5.34e6, delta=1.06e6, s=2,
        ignore_unconstrained=False):
    """
    Returns the cumulative number of craters at least T years old on Venus that
    were caused by an impactor of mass m.  Assumes that all impactors forms
    craters.  Alpha, beta, gamma, and delta, are parameters rescaled from
    Robbins 2013 that describe the cratering on Venus.  S is the power law
    parameter for mass scaling If ignore_constrained is false or not set, will
    throw an error if t any t is older than 3.92 Ga

    Inputs
    t:                          Gyr
    ignore_unconstrained:       Boolean
    alpha, beta, gamma, delta   ? (float) (Look at robbins to figure out units)
    s                           dimentionless

    Output
    Cumulative Craters          km^-2

    Source: Robbins (2013)
    """
    if (np.any(T > 3.92) and not ignore_unconstrained): raise

    return ((alpha * (np.exp(beta * T) - 1) + gamma * T + delta * (T**2))
            * (m/5.47e7)**-0.59)

def N_m_total(T_start, T_end, m):
    """
    Returns the total number of impactors of at least mass m that impacted Venus
    between T_start and T_end Ga.

    Inputs
    T_start                     Gya
    T_end                       Gya
    m                           kg

    Output
    Cumulative Craters          km^-2
    """
    return N_m(T_start, m) - N_m(T_end, m)


def burst_height(scale_height, v_0, yield_stress, pancake_factor, theta, radius,
                 m_0, drag_coeff):
    """
    Returns the height of catastrophic fragmentation of a impactor incident on
    a planet.  Scale_height is the scale height of the planetary atmosphere,
    v_0 is the incident velocity, yield_stress is the yield_stress of the
    impactor, pancake_factor is the "pancake factor" as described in
    Robertson et al. 2019, theta is the incident angle, radius is the initial
    radius of the impactor, m_0 is the initial mass of the impactor, and
    drag_coeff is the drag coefficient of the impactor.  If the return value is
    negative, that indicates a catastrophic fragmentation event did not happen.

    Inputs:
    scale_height            m
    v_0                     m/s
    yield_stress            Pa
    pancake_factor          dimentionless
    theta                   dementionless (radians)
    radius                  m
    m_0                     kg
    drag_coeff              dimentionless

    Output:
    Burst Height            m

    Source Robertson et al. 2019
    """

    E = 1/2 * m_0 * v_0**2


    return (scale_height * np.log(1/2 * rho_0 * v_0**2 / yield_stress)
            - 2 * scale_height
            * np.log(1 + pancake_factor*np.sin(theta)/(2*scale_height)
                     * np.sqrt(6 * E/(2*drag_coeff*np.pi*yield_stress*radius))))

@numba.njit
def atmospheric_density(scale_height, z):
    """
    Returns the density of the atmosphere assuming simple exponential decay
    as a function of altitude, z.

    Inputs:
    scale_height            m
    z                       m

    Output
    atmospheric_density     kg/m^3

    Source:
    Standard Equation (Robertson et al. 2019)
    """

    return rho_0 * np.exp(-z/scale_height)




def atmospheric_drag(drag_coeff, rho_a, A, v, m, g, theta):
    """
    Returns the instantaneous acceleration due to atmospheric drag as a function
    of the drag coefficient, drag_coeff, the density of the atmosphere,
    rho_a, the cross_sectional area of the object, A, the
    velocity of the object, v, the mass of the object, m, the gravitational
    field strength, g, and the trajectory angle relative to the ground, theta.

    Inputs:
    drag_coeff              dimentionless
    atmospheric_density     kg/m^3
    A                       m^2
    v                       m/s
    m                       kg
    g                       m/s^2
    theta                   dimentionless (radians)

    Outputs:
    atmospheric drag:       m/s^2

    Source:
    Melosh 1989
    """

    return -drag_coeff * rho_a * A * v**2 / m + g * np.sin(theta)

@numba.njit
def atmospheric_effects(t, X, m, drag_coeff, rho_imp, lift_coeff, planet_id):
    """
    Takes in a state vector X = [h, v, theta] and returns the instantaneos rate
    of change in altitude, h, velocity, v, and angle of trajectory relative to
    ground, theta, as a function of the mass of the object, m, the drag
    coefficent, drag_coeff, the density of the impactor, rho_imp, the lift
    coefficient, lift_coeff, and the planet id (Earth = 3, Venus = 2 etc.) t is
    a dummy variable for the integrator methods.

    Input:
    t               does not matter
    X               vector ([h, v, theta])
        h           m
        v           m/s
        theta       dimentionless (radians)
    m               kg
    drag_coeff      dimentionless
    rho_imp         kg/m^3
    lift_coeff      dimentionless
    planet_id       dimentionless

    Output:
    dhdt            m/s
    dvdt            m/s^2
    dthetadt        /s    (radians/sec)

    """


    h, v, theta = X
    g = GRAV_CONSTANT*VENUS_MASS/((VENUS_RADIUS + h)**2)
    r = (3 * m / (4 * rho_imp * np.pi))**(1/3)
    A = np.pi * r**2

    atm_rho = atmospheric_density(8000, h)


    dhdt     = -v*np.sin(theta)
    dvdt     = -(drag_coeff*atm_rho*A*(v**2)/m + g * np.sin(theta))



    dthetadt = ((g * np.cos(theta) / v)- (lift_coeff * atm_rho * A * v/(2*m))
                - (v * np.cos(theta)/(VENUS_RADIUS + h)))

    if h <= 0: return (0, 0, 0)
    if dhdt > 0: dhdt = 0
    if dhdt == 0: dvdt = 0
    return (dhdt, dvdt, dthetadt)

