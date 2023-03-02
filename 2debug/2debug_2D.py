#!/usr/bin/env python3

import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
import argparse
import pathlib
import h5py as h5

FRAME_X = -1.5
FRAME_Y = FRAME_X
HEIGHT = np.abs(2*FRAME_X)
WIDTH = HEIGHT

UPPER_LIM_X = 2.0
LOWER_LIM_X = -UPPER_LIM_X
UPPER_LIM_Y = UPPER_LIM_X
LOWER_LIM_Y = -UPPER_LIM_X

MARKER_SIZE_PAR = 100
MARKER_SIZE_NN = 50

def createPossitionPlot(h5File, outDir, flagGrad, flagVel, iNNL):
    data = h5.File(h5File, 'r')
    pos = data["x"][:]
    rho = data["rho"][()]
    data.close()

    fig, ax = plt.subplots()
    rhoPlt = ax.scatter(pos[:,0], pos[:,1], c=rho, s=MARKER_SIZE_PAR)

    # marks NNL for particle i
    if iNNL > -1 and "Ghosts" not in str(h5File):
        marksNN(h5File, iNNL, pos, ax)

    # set domain limits
    if "Ghosts" not in str(h5File):
        ax.set_xlim((LOWER_LIM_X, UPPER_LIM_X))
        ax.set_ylim((LOWER_LIM_Y, UPPER_LIM_Y))

    # Create a Rectangle patch
    rect = patches.Rectangle((FRAME_X, FRAME_Y), HEIGHT , WIDTH, alpha=0.1, facecolor="black")

    # Add the patch to the Axes
    ax.add_patch(rect)

    fig.colorbar(rhoPlt, ax=ax)
    
    plt.title(r"Color coded density $\rho$")
    plt.xlabel("$x$")
    plt.ylabel("$y$")
    plt.tight_layout()
    plt.grid()
    print("\tSaving figure to", outDir + "/" + pathlib.Path(h5File).stem + ".png")
    plt.savefig(outDir + "/" + pathlib.Path(h5File).stem + ".png")
    plt.close()

def marksNN(h5File, iNNL, pos, ax):
    data = h5.File(str(h5File).replace(".h5", "NNL.h5"), 'r')
    posNNL = data["nnlPrtcls"+str(iNNL)][:]
    data.close()

    ax.scatter(pos[iNNL,0], pos[iNNL,1], s=MARKER_SIZE_NN, marker='x', color='b')
    ax.scatter(posNNL[:,0], posNNL[:,1], s=MARKER_SIZE_NN, marker='x', color='r')

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
        
