"""
sim.py
author: Wolf Cukier
Simulates the collisional timescale of meteroid ejecta from Venus
"""
from timescale.constants import *
from timescale.particles import initial_state
from timescale.collisions import collision_probability
import rebound
import numpy as np
import sys
import json

with open("data/planets.json") as f:
    planets = json.load(f)
with open("config.json") as f:
    config = json.load(f)
    
N_ACTIVE = len(config["planets"]) + 1
YEAR_STEP = config["YEAR_STEP"]


def initialize(max_years, n):
    """
    Initialize the rebound simulation which will later be run for max_years.
    Returns the simulation and an array to log the output in.  Prints success to
    stderr.
    
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
    print("Simulation Initialized", flush=True, file = sys.stderr)
    
    return sim, np.zeros((int(max_years/YEAR_STEP)+1,n,4))
    

def add_particles(sim, n, v_inf, start = 0, end = -1, states = -1):
    """
    Add the n particles, along with the Sun, Venus, Earth, and Jupiter, to the
    simulation.  The particles are added after the major solar system bodies.
    If start and end are specified, only particles with array bounds between 
    start and end are actually added.  If states are not provided adds particles
    with a velocity of v_inf, otherwise uses the provided states. Prints sucess
    to stderr.
    
    Units
    sim:           (Rebound Simulation)
    n:             dimentionless
    v_inf:         m/s
    start:         dimentionless
    end:           dimentionless
    return:        void
    """
    if (end < 0): end = n
    
    if (np.all(states == -1)): 
        states = initial_state(n, v_inf, planets = config["planets"], 
                               source_id=config["source"])

    sim.add(m = MASS_SUN)
    sim.move_to_hel()
    
    for i , p in enumerate(config["planets"]):
        sim.add(m = planets[p]["mass"], 
                x = states[i,0], y = states[i,1], z = states[i,2],
                vx = states[i,3], vy = states[i,4], vz = states[i,5])
    
    for i in range(3 + start, end+3):
        sim.add(x = states[i,0], y = states[i,1], z = states[i,2],
                vx = states[i,3], vy = states[i,4], vz = states[i,5],
                hash = f"{i-3-start}")
        
    sim.move_to_com()
    sim.n_active = N_ACTIVE

    print(f"Planets and and {n} particles have been added.",
         flush = True, file=sys.stderr)
    
    return
    
def step(sim, year):
    """
    Advances the simulation sim forward to year years after simulation start. 
    Prints the current execution year to stderr.
    
    Units
    sim:        (Rebound Simulation)
    years:      years
    return:     void
    """
    sim.integrate(year * SEC_PER_YEAR, exact_finish_time = 0)
    print(f"=====Year {year}=====", flush = True, file = sys.stderr)
    return

    
def remove_particles(sim, n_removed):
    """
    Removes all particles that are either sun-grazing (<.01 au) or escaping the
    system (>10 au).  Prints the number of removed particles, n_removed, to 
    stderr
    
    Units
    sim:        (Rebound Simulation)
    n_removed:  dimentionless
    return:     dimentionless
    """
    
    sim.integrator_synchronize()
    for p in sim.particles[N_ACTIVE:]:
        h = p.hash
        d = sim.particles[0] ** p
        
        if d < .01*AU:
            sim.remove(hash = h)
            n_removed += 1
            print(f"Removed sun-grazing particle. \
                  {n_removed} particles removed", flush = True, 
                  file = sys.stderr)
            sim.ri_whfast.recalculate_coordinates_this_timestep = 1
            
        elif d > 10*AU:
            sim.remove(hash = h)
            n_removed += 1
            print(f"Removed escaping particle. \
                  {n_removed} particles removed", flush = True,
                  file = sys.stderr)
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
            prob = collision_probability(o.a, o.e, o.inc, 
                                         YEAR_STEP*SEC_PER_YEAR,
                                         planets[config["target"]])

            logger[step, h, :] = [o.a, o.e, o.inc, prob]
        except Exception as e:
            print(e)
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
    
    
def simulate(n, max_years, v_inf):
    """
    Simulates n particles launched 300 planet radii from the source planet at a 
    velocity of v_inf for max_years amount of time.  Returns a log of the 
    semi-major axis, eccentricity, inclination, and probability of hitting earth
    for each particle as a numpy array. If start and end are specified, only 
    particles with array bounds between start and end are actually added.
    
    Units
    n:           dimentionless
    max_years:   years
    v_inf:       m/s
    start:       dimentionless
    end:         dimentionless
    return:      m, dimentionless, degrees, dimentionless (max_years+1, n, 4)
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
        
    return logger


def sim_set_states(n, max_years, v_inf, start, end, states):
    """
    Simulates n particles launched 300 planet radii from the source planet at a 
    velocity of v_inf for max_years amount of time.  Returns a log of the 
    semi-major axis, eccentricity, inclination, and probability of hitting earth
    for each particle as a numpy array. If start and end are specified, only 
    particles with array bounds between start and end are actually added.
    Simulation uses the states provided but only the start-th to end-th 
    non-active particles
    
    Units
    n:           dimentionless
    max_years:   years
    v_inf:       m/s
    start:       dimentionless
    end:         dimentionless
    states:      [:, 0:3] m, [:, 3:] m/s
    return:      m, dimentionless, degrees, dimentionless (max_years+1, n, 4)
    """
    if (end < 0): end = n
    
    sim, logger = initialize(max_years, end-start)
    add_particles(sim, n, v_inf, start=start, end=end, states=states)
    
    year = 0
    n_removed = 0
    while(year <= max_years):
        step(sim, year)
        n_removed += remove_particles(sim, n_removed)
        log(sim, logger, end-start, year)
        year += YEAR_STEP
        
    return logger


if (__name__ == "__main__"):
    logger = simulate(int(sys.argv[1]), int(sys.argv[2]), float(sys.argv[3]))
    write_log(logger, float(sys.argv[3]), int(sys.argv[4]))
