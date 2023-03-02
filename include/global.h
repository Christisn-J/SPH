//
// Created by Christian Jetter on 23.01.23
//

#ifndef GLOBAL_H
#define GLOBAL_H

#include <string>
#include <namespace.h>
#include <parameter.h>

#define DISABLE -1

#define MAX_GHOSTS_PER_PARTICLE int(2*DIM-1)

struct Configuration {
    std::string initFile;
    std::string outDir;
    double timeStep;
    double timeStart;
    double timeEnd;
    int storeFrequency;
    int maxInteractions;
    double h; // smoothing lenght
    double gamma; // adiabatic index

#if BOUNDARIES != TRANSPARENT
    int maxGhostInteractions;
    double boxLimits[2 * DIM];
#endif // NOT TRANSPARENT 
};

#endif // GLOBAL_H
