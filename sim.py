"""
sim.py
author: Wolf Cukier
Simulates the collisional timescale of meteroid ejecta from Venus
"""
from constants import *
from particles import initial_state
from collisions import collision_probability
import rebound
import numpy as np
import sys


N_ACTIVE = 4
YEAR_STEP = 5000


def initialize(max_years, n):
    """
    Initialize the rebound simulation which will later be run for max_years.
    Returns the simulation and an array to log the output in
    
    Units
    return:        (Rebound Simulation), dimentionless (max_years+1, n, 4)
    """
    sim = rebound.Simulation()
    sim.units = ("s", "m", "kg")
    sim.integrator = "mercurius"
    sim.ri_whfast.safe_mode = 0
    sim.ri_whfast.corrector = 11
    sim.dt = 1e4
    sim.ri_ias15.min_dt = 1e-4 * sim.dt
    sim.testparticle_type = 0
    
    
    return sim, np.zeros((int(max_years/YEAR_STEP)+1,n,4))
    

def add_particles(sim, n, v_inf):
    """
    Add the n particles, along with the Sun, Venus, Earth, and Jupiter, to the
    simulation.  The particles are added after the major solar system bodies
    
    Units
    sim:           (Rebound Simulation)
    n:             dimentionless
    v_inf:         m/s
    return:        void
    """
    states = initial_state(n, v_inf, planets = ["2", "3", "5"])
    
    sim.add(m = MASS_SUN)
    sim.move_to_hel()
    
    sim.add(m = MASS_VENUS, x = states[0,0], y = states[0,1], z = states[0,2], 
            vx = states[0,3], vy = states[0,4], vz = states[0,5])
    sim.add(m = MASS_EARTH, x = states[1,0], y = states[1,1], z = states[1,2], 
            vx = states[1,3], vy = states[1,4], vz = states[1,5])
    sim.add(m = MASS_VENUS, x = states[2,0], y = states[2,1], z = states[2,2], 
            vx = states[2,3], vy = states[2,4], vz = states[2,5])
    
    for i in range(3, n+3):
        sim.add(x = states[i,0], y = states[i,1], z = states[i,2],
                vx = states[i,3], vy = states[i,4], vz = states[i,5],
                hash = f"{i-3}")
        
    sim.move_to_com()
    sim.n_active = 4
    
    return
    
def step(sim, year):
    """
    Advances the simulation forward to year years after simulation start
    
    Units
    sim:        (Rebound Simulation)
    years:      years
    return:     void
    """
    sim.integrate(year * SEC_PER_YEAR, exact_finish_time = 0)
    print(f"=====Year {year}=====")
    return

    
def remove_particles(sim, n_removed):
    """
    Removes all particles that are either sun-grazing (<.01 au) or escaping the
    system (>10 au).
    
    Units
    sim:        (Rebound Simulation)
    n_removed:  dimentionless
    """
    
    sim.integrator_synchronize()
    for p in sim.particles[N_ACTIVE:]:
        h = p.hash
        d = sim.particles[0] ** p
        
        if d < .01*AU:
            sim.remove(hash = h)
            n_removed += 1
            print(f"Removed sun-grazing particle. \
                  {n_removed} particles removed")
            sim.ri_whfast.recalculate_coordinates_this_timestep = 1
            
        elif d > 10*AU:
            sim.remove(hash = h)
            n_removed += 1
            print(f"Removed escaping particle. \
                  {n_removed} particles removed")
            sim.ri_whfast.recalculate_coordinates_this_timestep = 1
            
    return n_removed
        
def log(sim, logger, n, year):
    """
    Stores the semi-major axis, eccentricity, and inclination for each of n
    massless particles in sim into log, in the corresponding location based on
    year.
    
    Units
    sim           (Rebound Simulation)
    log           m, dimentionless, degrees, dimentionless (max_years+1, n, 4)
    n             dimentionless
    year          years
    return        void
    """
    
    
    step = int(year/YEAR_STEP)
    for h in range(n):
        try:
            p = sim.particles[f"{h}"]
            o = p.calculate_orbit(primary = sim.particles[0])
            prob = collision_probability(o.a, o.e, o.inc, YEAR_STEP)
            logger[step, h, :] = [o.a, o.e, o.inc, prob]
        except:
            logger[step, h, :] = [np.nan, np.nan, np.nan, np.nan]

    return

def write_log (logger, v_inf, run_num):
    """
    Writes the data contained in logger to a .npy file on disk.  File name is
    f"v{int(v_inf/1000)}-r{run_num}".
    
    Units are irrelevant
    returns void
    """
    np.save(f"output/v{int(v_inf/1000)}-r{run_num}", logger)
    return
    
    
def simulate(n, max_years, v_inf, run_num):
    """
    Simulates n particles launched 300 Venus radii from Venus at a velocity of
    v_inf for max_years amount of time.  Logs the semi-major axis, eccentricity,
    inclination, and chance of hitting earth for each particle in 
    f"{v_inf}_{run_num}.npy".
    
    Units
    n:           dimentionless
    max_years:   years
    v_inf:       m/s
    run_num:     dimentionless
    return:      void
    """
    
    sim, logger = initialize(max_years, n)
    add_particles(sim, n, v_inf)
    
    year = 0
    n_removed = 0
    while(year <= max_years):
        step(sim, year)
        n_removed += remove_particles(sim, n_removed)
        log(sim, logger, n, year)
        year += YEAR_STEP
    
    
    write_log(logger, v_inf, run_num)
    return

if (__name__ == "__main__"):
    simulate(int(sys.argv[1]), int(sys.argv[2]), float(sys.argv[3]), 
             int(sys.argv[4]))
    
    
    
    
    
