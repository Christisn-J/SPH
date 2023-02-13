#include "../include/Particles.h"

Particles::Particles(int numParticles) : N { numParticles }{
    // allocate memory
    m = new double[numParticles];
    u = new double[numParticles];

    rho = new double[numParticles];
    rhoGrad = new double[numParticles][DIM];

    P = new double[numParticles];
    PGrad = new double[numParticles][DIM];

    x = new double[numParticles];
    vx = new double[numParticles];
    vxGrad = new double[numParticles][DIM];
    
    #if DIM == 2
        y = new double[numParticles];
        vy = new double[numParticles];
        vyGrad = new double[numParticles][DIM];
    #endif // 2D

    #if DIM == 3
        z = new double[numParticles];
        vz = new double[numParticles];
        vzGrad = new double[numParticles][DIM];
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

    #if DIM == 2
        delete[] y;
        delete[] vy;
        delete[] vyGrad;
    #endif // 2D

    #if DIM == 3
            delete[] z;
            delete[] vz;
            delete[] vzGrad;
    #endif // 3D
}

