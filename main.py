"""
main.py
author: Wolf Cukier
Simulates the timescale for meteor ejecta traveling from Venus to Earth. 
Multithreaded program, designed specifically to be run using slurm
"""

from timescale.sim import sim_set_states, write_log
from multiprocessing import Pool
from timescale.particles import initial_state
import os
import sys
import numpy as np
import json

if __name__ == "__main__":
    with open("config.json") as f:
        config = json.load(f)
    num_cores = int(os.getenv('SLURM_CPUS_PER_TASK'))
    PARTICLES_PER_RUN = config["PARTICLES_PER_RUN"]
    planets = config["planets"]
    source = config["source"]

    n = config["n"]
    max_years = config["max_years"]
    v_inf = config["v_inf"]
    try:
        run_num = int(sys.argv[1])
    except:
        print("No run number provided.  Defaulting to 0.", file=sys.stderr)
    
  
    param_array = []
    start = 0
    end = PARTICLES_PER_RUN
    states = initial_state(n, v_inf, planets = planets, source_id=source)
    while(end < n):
        param_array.append([n, max_years, v_inf, start, end, states])
        start += PARTICLES_PER_RUN
        end += PARTICLES_PER_RUN
        
    if(start < n):
        param_array.append([n, max_years, v_inf, start, n, states])
        
    with Pool(num_cores) as pool:
        res = pool.starmap(sim_set_states, param_array)
    
    write_log(np.array(res), v_inf, run_num)
    
