//
// Created by Johannes Martin on 21.09.22.
//

#ifndef PARTICLES_H
#define PARTICLES_H

#include <cmath>
#include "global.h"
#include "kernel.h"

class Particles {
public:
    Particles(int numParticles);
    ~Particles();
        
    int N;
    double *m, *u, *x, *vx, *rho, *P;
    double (*rhoGrad)[DIM], (*vxGrad)[DIM], (*PGrad)[DIM];

#if DIM >= 2
    double *y, *vy;
    double (*vyGrad)[DIM];
#endif // 2D

#if DIM == 3
    double *z, *vz;
    double (*vzGrad)[DIM];
#endif // 3D

private:
    bool ghosts; // tells if particl is virtuel or real

    int *nnl; // nearest neighbor list
    int *noi;  // number of interactionsdouble 
    
    double (*kernel)(const double&, const double&){ &Kernel::cubicSpline };
};

#endif //PARTICLES_H