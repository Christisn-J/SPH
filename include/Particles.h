//
// Created by Johannes Martin on 21.09.22.
//

#ifndef PARTICLES_H
#define PARTICLES_H

#include <cmath>
#include <limits>
#include <highfive/H5File.hpp>

#include "global.h"
#include "Logger.h"
#include "kernel.h"
#include "Domain.h"

class Particles {
public:
    Particles(int numParticles, Configuration config);
    ~Particles();
        
    int N, MAX_INTERACTIONS;
    int *cell; // cell in which particle at index resides
    double *m, *u, *rho, *P;
    double (*rhoGrad)[DIM], (*PGrad)[DIM];

    double *x, *vx, *vxDelta;
    double (*vxGrad)[DIM];
#if DIM >= 2
    double *y, *vy, *vyDelta;
    double (*vyGrad)[DIM];
#endif // 2D
#if DIM == 3
    double *z, *vz, *vzDelta;
    double (*vzGrad)[DIM];
#endif // 3D

    void assignParticlesAndCells(Domain &domain);
    void compNN(Domain &domain, const double &h);
    void compDensity(const double &h);
    void compPressure(const double &gamma);
    void gradient(double *f, double (*grad)[DIM]);

    double sumMass();
    double sumMomentum(const int r);
    double sumEnergy();
    
    void accelerate(const double &h);
    void damping(const double &h);

    void integrate(const double &t, const Domain &domain);

    void getDomainLimits(double *domainLimits);
    void save(std::string filename, double simTime);

private:
    bool ghosts; // tells if particl is virtuel or real

    int *nnl; // nearest neighbor list
    int *noi; // number of interactionsdouble 
    
    double (*W)(const double&, const double&){ &KernelR::cubicSpline };
    double (*dW)(const double&, const double&){ &NablaKernelR::cubicSpline };

    int getNNLidx(int i, int n);

    double distance(const int i, const int n);

    double sumMomentumX();
#if DIM >= 2
    double sumMomentumY();
#endif // 2D
#if DIM == 3
    double sumMomentumZ();
#endif // 3D
};

#endif //PARTICLES_H