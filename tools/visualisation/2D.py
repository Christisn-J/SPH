#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import argparse
import pathlib
import h5py

def createPossitionPlot(h5File, outDir, flagGrad, flagVel, flagNN):
    data = h5py.File(h5File, 'r')
    pos = data["x"][:]
    rho = data["rho"][()]
    data.close()

    fig, ax = plt.subplots()
    rhoPlt = ax.scatter(pos[:,0], pos[:,1], c=rho, s=100.)

    ax.set_xlim((-5, 5))
    ax.set_ylim((-5, 5))

    fig.colorbar(rhoPlt, ax=ax)

    plt.xlabel("$x$")
    plt.ylabel("$y$")
    plt.tight_layout()
    plt.grid()
    print("\tSaving figure to", outDir + "/" + pathlib.Path(h5File).stem + ".png")
    plt.savefig(outDir + "/" + pathlib.Path(h5File).stem + ".png")
    plt.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot density of results from Kelvin-Helmholtz test case.")
    parser.add_argument("--inDir", "-d", metavar="string", type=str, help="output directory of simulation", required=True)
    parser.add_argument("--outDir", "-o", metavar="string", type=str, help="output directory for generated plots", default="out")
    parser.add_argument("--Gradient", "-g", action="store_true", help="plot density gradients")
    parser.add_argument("--Ghosts", "-G", action="store_true", help="also plot ghost cells in an extra file")
    parser.add_argument("--pressure", "-P", action="store_true", help="plot pressure instead of density")
    parser.add_argument("--energy", "-u", action="store_true", help="plot internal energy instead of density")
    parser.add_argument("--Velocity", "-v", action="store_true", help="plot velocity")
    parser.add_argument("--NN", "-i", metavar="int", type=int, help="plot NNL for particles i", default=-1)

    args = parser.parse_args()
    
    print("Examining files in", args.inDir, "...")
    
    for h5File in pathlib.Path(args.inDir).glob('*.h5'):
        if "NNL" not in str(h5File):
            if (args.Ghosts) or (not "Ghost" in str(h5File)): 
                print("\t", h5File)
                createPossitionPlot(h5File, args.outDir, args.Gradient, args.Velocity, args.NN)
        
