#ifndef PARTICLES_H
#define PARTICLES_H

#include "global.h"

class Particles {
    public:
        Particles(int numParticles);
        ~Particles();

        int N;
        double *m, *u, *x, *vx, *rho, *P;
        double (*rhoGrad)[DIM], (*vxGrad)[DIM], (*PGrad)[DIM];

        #if DIM == 2
            double *y, *vy;
            double (*vyGrad)[DIM];
        #endif // 2D

        #if DIM == 3
            double *z, *vz;
            double (*vzGrad)[DIM];
        #endif // 3D
};
#endif //PARTICLES_H