"""
main.py
author: Wolf Cukier
Simulates the timescale for meteor ejecta traveling from Venus to Earth. 
Multithreaded program, designed specifically to be run using slurm
"""

from sim import simulate, write_log
from multiprocessing import Pool
from particles import initial_state
import os
import sys
import numpy as np

if __name__ == "__main__":
    num_cores = int(os.getenv('SLURM_CPUS_PER_TASK'))
    PARTICLES_PER_RUN = 10
    
    
    n = int(sys.argv[1])
    max_years = int(sys.argv[2])
    v_inf = float(sys.argv[3])
    run_num = int(sys.argv[4])
    
    
    param_array = []
    start = 0
    end = PARTICLES_PER_RUN
    states = initial_state(n, v_inf, planets = ["2", "3", "5"])
    while(end < n):
        param_array.append([n, max_years, v_inf, start, end, states])
        start += PARTICLES_PER_RUN
        end += PARTICLES_PER_RUN
        
    if(start < n):
        param_array.append([n, max_years, v_inf, start, n, states])
        
    with Pool(num_cores) as pool:
        res = pool.starmap(sim_set_states, param_array)
    
    write_log(np.array(res), v_inf, run_num)

    
