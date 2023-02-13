#!/usr/bin/env python3

import numpy as np
from itertools import permutations
import h5py

DIM = 2
N = 10

def cube(center, length):
    pos = np.zeros((N, DIM))
    vel = np.zeros((N, DIM))
    mass = np.zeros(N)

    x = np.linspace(center, length, N)
    perm = np.array(permutations(x, DIM))

    for n in range(N):
        mass[n] = 1
        for i in range(DIM):
            vel[n,i] = float(0)
            pos[n,i] = x[n]

    return pos, vel, mass
            


def main(): 
    f= cube
    
    h5f = h5py.File("./N{}.h5".format(N), "w")
    pos, vel, mass = f(0,1)

    print("Writing to HDF5 file ...")
    h5f.create_dataset("x", data=pos)
    h5f.create_dataset("v", data=vel)
    h5f.create_dataset("m", data=mass)

    h5f.close()

    print("Finished!")

if __name__ == "__main__":
    main()