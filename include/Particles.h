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

enum typeNN {
    PROTFORCE,
};

class Particles {
public:
    Particles(int numParticles, Configuration config,  bool ghosts = false);
    ~Particles();
        
    int N, MAX_INTERACTIONS;
    double *m, *u, *x, *vx, *vxDelta, *rho, *P;
    double (*rhoGrad)[DIM], (*vxGrad)[DIM], (*PGrad)[DIM];

#if DIM >= 2
    double *y, *vy, *vyDelta;
    double (*vyGrad)[DIM];
#endif // 2D

#if DIM == 3
    double *z, *vz, *vzDelta;
    double (*vzGrad)[DIM];
#endif // 3D

    void compNN(const double &h);
    void compDensity(const double &h);
    void compPressure(const double &gamma);
    void gradient(double *f, double (*grad)[DIM]);

    double sumMass();
    double sumMomentum(const int r);
    double sumEnergy();
    
    void accelerate(const double &h);
    void damping(const double &h);

    void integrate(const double &t);

    void getDomainLimits(double *domainLimits);
    void checkBoundary(const Domain &domain);
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