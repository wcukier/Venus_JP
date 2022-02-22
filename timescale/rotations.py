"""
rotations.py
author: Wolf Cukier

Some helper functions that involve matricies
"""
import numpy as np

def R_m(L):
    """
    Takes a 3D vector, L, and generates a rotation matrix, M, such that
    ML = e_3|L| where e_3 is zhat and |L| is the norm of L.  Generates the axis
    angle representation of the transformation and then transforms it to matrix
    form.
    """
    L_hat = L / np.linalg.norm(L)
    e3 = np.array([0, 0, 1])
    e = np.cross(L_hat, e3)
    u = e / np.linalg.norm(e)
    theta = np.arccos(np.dot(L_hat, e3))
    return np.array([
        [
            np.cos(theta) + u[0]**2*(1-np.cos(theta)), 
            u[0]*u[1]*(1- np.cos(theta)) - u[2]*np.sin(theta), 
            u[0]*u[2]*(1-np.cos(theta)) + u[1]*np.sin(theta)
         ],
        [
            u[1]*u[0]*(1-np.cos(theta)) + u[2]*np.sin(theta),
            np.cos(theta) + u[1]**2*(1-np.cos(theta)),
            u[1]*u[2]*(1-np.cos(theta)) - u[0]*np.sin(theta)
        ],
        [
            u[2]*u[0]*(1 - np.cos(theta)) - u[1]*np.sin(theta), 
            u[2]*u[1]*(1- np.cos(theta)) + u[0]*np.sin(theta),
            np.cos(theta) + u[2]**2*(1-np.cos(theta))
        ]
    ])