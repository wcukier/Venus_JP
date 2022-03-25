"""
spalling.py
author: Wolf Cukier

A set of equations that describe the number and size of fragments launched into
space as a result of an impact event
"""

import numpy as np

def mass_ej(P_max, rho_T, C_T, v_imp, v_0, m_imp):
    """
    Returns the mass of unshocked particles (subject to a pressure less than
    P_max) ejected at a speed of at less than V_0.  Rho_T is the density of the
    target (the planet), C_T is the longitudinal speed of sound in the target,
    v_imp is the impactor velocity normal to the surface of the target at time
    of impact, m_imp is the mass of the impactor at the time of impact.

    Inputs:
    P_max           Pa
    rho_T           kg/m^3
    C_T             m/s
    v_imp           m/s
    v_0             m/s
    m_imp           kg

    Output:
    m_ej            kg

    Source: Armstrong et al. 2016 and Melosh 1985
    """

    return (0.75 * (P_max /(rho_T * C_T * v_imp))
            * ((v_imp/(2*v_0))**(5/3) - 1)
            * m_imp )

def mean_size(rho_imp, T, rho_T, v_ej, v_imp, m_imp):
    """
    Returns the mean particle size (diameter?) launched in a spallation event as
    a function of the impactor density, rho_imp, the dynamic tensile strength of
    the target, T, the target density, rho_T, the limiting ejection velocity,
    v_ej, the impactor velocity, v_imp, and the impactor mass, m_imp.

    Inputs:
    rho_imp         kg/m^3
    T               Pa
    rho_T           kg/m^3
    v_ej            m/s
    v_imp           m/s
    m_imp           kg

    Output:
    size            m

    Source: Grady & Kipp 1980, Melosh 1985
    """

    return ((6  * (T**3))
            / (np.pi * rho_imp * (rho_T**3)
               * (v_ej**2) * (v_imp **4)) * m_imp)**(1/3)


def max_size(m_w, mean_s):
    """
    Returns the maximum fragment size launched in a spalling event as a function
    of the Weibulls constant, m_w, and the mean size of the launched particles,
    mean_size

    Inputs:
    m_w             dimentionless
    mean_size       m

    Output:
    max_size        m

    Source
    Melosh et al. 1992
    """

    return (m_w + 2) / 3 * mean_s

def n_launched(m_ej, mean_size, rho_T):
    """
    Returns the total number of particles launched by the spalling process as a
    function of the ejected mass, m_ej, the mean particle size, mean_size, and
    the target density, rho_T

    Inputs:
    m_ej            kg
    mean_size       m
    rho_T           kg/m^3

    Outputs:
    n_launched       dimentionless

    Source: Mileikowski et al. 2000
    """

    return m_ej / (np.pi/6 * mean_size**3 * rho_T)

def size_distribution(limit_size, max_size, m_w):
    """
    Returns the cumulative fraction of spalled fragments with a size larger than
    limit_size as a function of the Weibull's constant, m_w, and the
    maximum size of the spalled materian.

    Inputs:
    limit_size      m
    max_size        m
    m_w             dimentionless

    Output:
    size_fraction   dimentionless

    Source:
    Melosh et al. 1992
    """

    return ((1 - (limit_size/max_size))**(m_w)
            * (1
               + m_w*(limit_size/max_size)
               + m_w*(m_w + 1)/2 * (limit_size/max_size)**2
               + m_w*(m_w+1)*(m_w+2)/6*(limit_size/max_size)**3))


def v_inf(v_ej, v_esc):
    return np.sqrt(v_ej**2 - v_esc**2)

# def atmospheric_escape_vel(v_esc, m_air, beta, C_T, P_max, rho_imp, m_imp):
#     return (v_esc
#             / (1
#                - ((m_air * beta * C_T * v_esc * (rho_imp * 4 * np.pi)**(1/3))
#                   / (2 * P_max * (3 * m_imp)**(1/3)))))


# def v_inf(v_ej, v_esc, m_air, beta, C_T, P_max, rho_imp, m_imp):
#     v_lim = atmospheric_escape_vel(v_esc, m_air, beta,
#                                    C_T, P_max, rho_imp, m_imp)
#     return np.sqrt(v_ej **2 - v_lim**2)

