#!/usr/bin/env python3

import argparse
import numpy as np
import h5py

# default values ---------------------------------------------------------------------------------------------------------------------
DIM = 2
ORIGIN = 0.0
LENGTH = 1.0

# functions ---------------------------------------------------------------------------------------------------------------------
def grid(center, length, n, dim, flag =False):
    pos = np.zeros((n, dim))
    vel = np.zeros((n, dim))
    mass = np.zeros(n)

    
    for i in range(n):
        mass[i] = 1
        for j in range(dim):
            vel[i,j] = float(1)

    i = 0
    for x in np.linspace(center-length,length,int(np.sqrt(n))):
        if dim == 1:
            if flag:
                pos[i] = [x] 
            else:
                pos[i] = 2*length*np.random.random(1) + center-length
            i += 1
        else:
            for y in np.linspace(center-length,length,int(np.sqrt(n))):
                if dim == 2:
                    if flag:
                        pos[i] = [x,y] 
                    else:
                        pos[i] = 2*length*np.random.random(2) + center-length
                    i += 1
                else:
                    for z in np.linspace(center-length,length,int(np.sqrt(n))):
                        if dim == 3:
                            if flag:
                                pos[i] = [x,y,z] 
                            else:
                                pos[i] = 2*length*np.random.random(3) + center-length
                            i += 1
                        else:
                            print(ERROR)
                   
    return pos, vel, mass

# main function ---------------------------------------------------------------------------------------------------------------------
def main(): 
    parser = argparse.ArgumentParser(description="Create an initial condition HDF5 file for a 2D Kelvin-Helmholtz test case.")
    parser.add_argument("--particles", "-N", metavar="int", type=int, help="number of particles", required=True)
    parser.add_argument("--origin", "-o", metavar="float", type=float, help="Origin of the Grid", default=ORIGIN)
    parser.add_argument("--length", "-l", metavar="float", type=float, help="half length of domain", default=LENGTH)
    parser.add_argument("--dim", "-D", metavar="float", type=int, help="Number of dimensions", default=DIM)
    parser.add_argument("--regularGrid", "-g", action="store_true")

    args = parser.parse_args()

    h5f = h5py.File("./2debug_N{}_{}D.h5".format(args.particles, args.dim), "w")
    if args.regularGrid:
        print("Generating regular grid positions")
    else:
        print("Generating random positions")
    pos, vel, mass = grid(args.origin, args.length, args.particles, args.dim, args.regularGrid)

    print("Writing to HDF5 file ...")
    h5f.create_dataset("x", data=pos)
    h5f.create_dataset("v", data=vel)
    h5f.create_dataset("m", data=mass)

# TODO: not ISOTHERM (Energy chancing)
    u = np.zeros(args.particles)
    h5f.create_dataset("u", data=u)

    h5f.close()

    print("Finished!")

#---------------------------------------------------------------------------------------------------------------------
if __name__ == "__main__":
    main()