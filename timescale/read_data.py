"""
read_data.py
author: Wolf Cukier
Reads and processes the data produced by the Venus_JP simulation.
"""
import numpy as np
import json
from tqdm import tqdm
from timescale.constants import *
from timescale.collisions import collision_probability

def load_data(v_inf, run_num, path="output/"):
    """
    Loads the output file specified by v_inf and run and returns the data as a
    numpy array.
    
    Units
    v_inf:      m/s
    run:        dimentionless
    """
    
    data = np.load(f"{path}v{int(v_inf/1000)}-r{run_num}.npy")
    return collate_data(data)


def collate_data(data):
    """
    Takes data, a 4d np array which is the output of the timescale simulation
    and returns the data in a 3d array such that the different cpu cores are now
    all along the same axis.

    Units
    data:       m, dimentionless, degrees, dimentionless
    """
    new_data = data[0]
    for i in range(1, data.shape[0]):
        new_data = np.hstack((new_data, data[i]))

    return new_data


def cume_prob(data):
    """
    Takes data, an np array which is the output of the timescale simulation, and
    returns a np array with the cumulative probability.

    Units
    data:       m, dimentionless, degrees, dimentionless
    """

    part_probs = data[:,:,3]
    part_probs[np.isnan(part_probs)] = 0
    cume_part_probs = np.ones((part_probs.shape[0]+1, part_probs.shape[1]))

    for t in range(part_probs.shape[0]):
        for p in range(part_probs.shape[1]):
            if (~np.isnan(part_probs[t][p]) and part_probs[t][p] > 0
                and part_probs[t][p] < 1):
                cume_part_probs[t+1][p] = cume_part_probs[t][p] * (
                    1-part_probs[t][p])
            else:
                cume_part_probs[t+1][p] = cume_part_probs[t][p]

    cume_probs = np.ones(cume_part_probs.shape[0])

    for i in range(cume_part_probs.shape[0]):
        cume_probs[i] = np.mean(cume_part_probs[i])

    return  1 - cume_probs

def recalculate_data(data, target_id, YEAR_STEP = 5000):
    """
    Takes in data and the key of the target, target_id, and returns the
    cume_prob of the data with that target

    Input
    data:         m, dimentionless, degrees, dimentionless
    target:       string "2"-"8"

    Output
    cume_probs:   array, dimentionless (probability)
    """


    with open("../data/planets.json") as f:
        planets = json.load(f)
        target = planets[target_id]

    for t in tqdm(range(data.shape[0])):
        for p in range(data.shape[1]):
            data[t,p,3] = collision_probability(data[t,p,0], data[t,p,1],
                                                data[t,p,2],
                                                YEAR_STEP*SEC_PER_YEAR,
                                                target)


    return cume_prob(data)

