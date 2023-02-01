#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import h5py

N=100
DIM = 2
PAR = N**DIM

def main():
    h5f = h5py.File("./out/N{}{}D.h5".format(N,DIM), 'r')
    pos = np.array(h5f.get('pos'))
    rho = np.array(h5f.get('rho'))
    h5f.close()

    x=np.zeros(PAR)
    y=np.zeros(PAR)
    val=np.zeros(PAR)
    for n in range(PAR):
        x[n] = pos[n][0]
        y[n] = pos[n][1]
        val[n] = rho[n]

    fig, ax = plt.subplots()
    im = ax.scatter(x, y,s=1, c=val)
    fig.colorbar(im, ax=ax, label=r"$\rho$")
    ax.grid()

    ax.set_xlabel("x")
    ax.set_ylabel("y")

    fig.savefig("./out/N{}{}D.png".format(N,DIM)) 

if __name__ == "__main__":
    main()
