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
// constructor ----------------------------------------------------------------------------------------------- 
    Particles(int numParticles, Configuration config,  bool ghosts = false);
    
// deconstructor --------------------------------------------------------------------------------------------- 
    ~Particles();

// object variablen ---------------------------------------------------------------------------------------------  
    int N, MAX_INTERACTIONS;
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

#if NNS == GRID
    int *cell; // cell in which particle at index resides
#endif // NNS GRID

#if BOUNDARIES != TRANSPARENT
// ghost variablen -------------------------------------------------------------------------------------------------
    int MAX_GHOST_INTERACTIONS;
#endif // NOT TRANSPARENT

// functions -------------------------------------------------------------------------------------------------
    void accelerate(const double &h);
    void damping(const double &h);
    void integrate(const double &t);
    void boundary(const Domain &domain);    
    void save(std::string filename, double simTime);

// NNS functions -------------------------------------------------------------------------------------------------------
void compNNS(Domain &domain, const double &h);
#if NNS == GRID
    void assignParticlesAndCells(Domain &domain);
#endif // NNS GRID

// boundary functions-------------------------------------------------------------------------------------------------------
    void getDomainLimits(double *domainLimits);

// compelation functions -------------------------------------------------------------------------------------------------------
    void compDensity(const double &h);
    void compPressure(const double &gamma);
    void gradient(double *f, double (*grad)[DIM]);

// sanity check functions ------------------------------------------------------------------------------------
    double sumMass();
    double sumMomentum(const int r);
    double sumEnergy();
    
#if BOUNDARIES != TRANSPARENT
// ghost functions -------------------------------------------------------------------------------------------------
    void createGhostParticles(Domain &domain, Particles &ghosts, const double &h);
    void compNNG(Domain &domain, const Particles &ghosts, const double &h);
    void compDensity(const Particles &ghosts, const double &h);
    void updateGhostState(Particles &ghosts);
    void saveNNL(std::string filename, const Particles &ghosts);
    void accelerate(const Particles &ghosts, const double &h);
    void damping(const Particles &ghosts, const double &h);
#endif // NOT TRANSPARENT
  
private:
// object variablen ---------------------------------------------------------------------------------------------  
    bool ghosts; // tells if particl is virtuel or real

    int *nnl; // nearest neighbor list
    int *noi; // number of interactionsdouble

#if BOUNDARIES != TRANSPARENT
// ghost variablen -------------------------------------------------------------------------------------------------
    int *nnlGhosts;
    int *noiGhosts;
    int *ghostMap;
    int *parent; // if class holds ghost particles the parent node is stored
#endif // NOT TRANSPARENT

// namespace ------------------------------------------------------------------------------------
    double (*W)(const double&, const double&){ &KernelR::cubicSpline };
    double (*dW)(const double&, const double&){ &NablaKernelR::cubicSpline };

// functions -------------------------------------------------------------------------------------------------
    double distance(const int i, const int n);

#if BOUNDARIES != TRANSPARENT
// ghost functions -------------------------------------------------------------------------------------------------
    double distance2G(const Particles &ghosts,const int i, const int ip);
#endif // NOT TRANSPARENT

 // helper functions ------------------------------------------------------------------------------------------  
    int getNNLidx(const int &i, const int &n);

    double sumMomentumX();
#if DIM >= 2
    double sumMomentumY();
#endif // 2D
#if DIM == 3
    double sumMomentumZ();
#endif // 3D
};

#endif //PARTICLES_H