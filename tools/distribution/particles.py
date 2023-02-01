#!/usr/bin/env python3

import numpy as np
import h5py

N=100
DIM = 2
PAR = N**DIM


def grid(center, length, randam =False):
    pos = np.zeros((PAR, DIM))
    vel = np.zeros((PAR, DIM))
    mass = np.zeros(PAR)

    
    for n in range(N):
        mass[n] = 1
        for i in range(DIM):
            vel[n,i] = float(0)

    n = 0
    shift = length/2
    for x in np.linspace(center-shift,shift,N):
            for y in np.linspace(center-shift,shift,N):
                #pos[n] = [x,y] # grid pos
                pos[n] = length*np.random.random(2) + center-shift # random pos
                n += 1

    return pos, vel, mass

def main(): 
    # parser arguments
    OFFSET = 0
    POINT = 0
    
    h5f = h5py.File("./out/N{}{}D.h5".format(N,DIM), "w")
    pos, vel, mass = grid(0,100)

    print("Writing to HDF5 file ...")
    h5f.create_dataset("pos", data=pos)
    h5f.create_dataset("v", data=vel)
    h5f.create_dataset("m", data=mass)
    rho = np.zeros(PAR)
    h5f.create_dataset("rho", data=rho)

    h5f.close()

    print("Finished!")

if __name__ == "__main__":
    main()