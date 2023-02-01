//
// Created by Christian Jetter on 23.01.23
//

#include "../include/Particles.h"

Particles::Particles(int nParticles) : N { nParticles }{
    // allocate memory
    m = new double[N];
    u = new double[N];

    rho = new double[N];
    rhoGrad = new double[N][DIM];

    P = new double[N];
    PGrad = new double[N][DIM];

    x = new double[N];
    vx = new double[N];
    vxGrad = new double[N][DIM];
    
#if DIM >= 2
    y = new double[N];
    vy = new double[N];
    vyGrad = new double[N][DIM];
#endif // 2D

#if DIM == 3
    z = new double[N];
    vz = new double[N];
    vzGrad = new double[N][DIM];
#endif // 3D
}

Particles::~Particles() {
    delete[] m;
    delete[] u;

    delete[] rho;
    delete[] rhoGrad;

    delete[] P;
    delete[] PGrad;

    delete[] x;
    delete[] vx;
    delete[] vxGrad;

#if DIM >= 2
    delete[] y;
    delete[] vy;
    delete[] vyGrad;
#endif // 2D

#if DIM == 3
        delete[] z;
        delete[] vz;
        delete[] vzGrad;
#endif // 3D

#if BOUNDARIES == PERIODIC
#endif
}

