#ifndef PARTICLES_H
#define PARTICLES_H

#include "global.h"

class Particles {
    public:
        Particles(int numParticles);
        ~Particles();

        int N;
        double *m, *u, *x, *y, *vx, *vy, *rho, *P;
        double (*rhoGrad)[DIM], (*vxGrad)[DIM], (*vyGrad)[DIM], (*vzGrad)[DIM], (*PGrad)[DIM];
        #if DIM == 3
            double *z, *vz;
        #endif // DIM
}
#endif //PARTICLES_H