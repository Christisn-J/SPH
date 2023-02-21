#!/usr/bin/env python3

import argparse
import numpy as np
import h5py

DIM = 2


def grid(center, length, N, flag =False):
    pos = np.zeros((N, DIM))
    vel = np.zeros((N, DIM))
    mass = np.zeros(N)

    
    for n in range(N):
        mass[n] = 1
        for i in range(DIM):
            vel[n,i] = float(0)

    i = 0
    for x in np.linspace(center-length,length,int(np.sqrt(N))):
            for y in np.linspace(center-length,length,int(np.sqrt(N))):
                if flag:
                    pos[i] = [x,y] 
                else:
                    pos[i] = 2*length*np.random.random(2) + center-length
                i += 1
    return pos, vel, mass

def main(): 
    parser = argparse.ArgumentParser(description="Create an initial condition HDF5 file for a 2D Kelvin-Helmholtz test case.")
    parser.add_argument("--particles", "-N", metavar="int", type=int, help="number of particles", required=True)
    parser.add_argument("--origin", "-O", metavar="float", type=float, help="Origin of the Grid", default=0.0)
    parser.add_argument("--domain", "-D", metavar="float", type=float, help="number of particles", default=1.0)
    parser.add_argument("--regularGrid", "-g", action="store_true")

    args = parser.parse_args()

    N = args.particles
    
    h5f = h5py.File("./kh_N{}_{}D.h5".format(N,DIM), "w")
    if args.regularGrid:
        print("Generating regular grid positions")
    else:
        print("Generating random positions")
    pos, vel, mass = grid(args.origin, args.domain, N, args.regularGrid)

    print("Writing to HDF5 file ...")
    h5f.create_dataset("x", data=pos)
    h5f.create_dataset("v", data=vel)
    h5f.create_dataset("m", data=mass)

    u = np.zeros(N)
    h5f.create_dataset("u", data=u)

    h5f.close()

    print("Finished!")

if __name__ == "__main__":
    main()