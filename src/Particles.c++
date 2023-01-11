#include "../include/Particles.h"

Particles::Particles(int numParticles) : N { numParticles }{
    // allocate memory
    m = new double[numParticles];
    u = new double[numParticles];
    rho = new double[numParticles];
    P = new double[numParticles];
    x = new double[numParticles];
    y = new double[numParticles];
    vx = new double[numParticles];
    vy = new double[numParticles];
    rhoGrad = new double[numParticles][DIM];
    vxGrad = new double[numParticles][DIM];
    vyGrad = new double[numParticles][DIM];
    vzGrad = new double[numParticles][DIM];
    PGrad = new double[numParticles][DIM];

    #if DIM == 3
        z = new double[numParticles];
        vz = new double[numParticles];
    #endif // DIM
}

Particles::~Particles() {
    delete[] m;
    delete[] u;
    delete[] P;
    delete[] rho;
    delete[] x;
    delete[] y;
    delete[] vx;
    delete[] vy;
    delete[] rhoGrad;
    delete[] vxGrad;
    delete[] vyGrad;
    delete[] vzGrad;
    delete[] PGrad;

    #if DIM == 3
            delete[] z;
            delete[] vz;
    #endif
}

