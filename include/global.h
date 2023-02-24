//
// Created by Christian Jetter on 23.01.23
//

#ifndef GLOBAL_H
#define GLOBAL_H

#include <parameter.h>
#include <Boundary.h>
#include <string>

/// TODO: remove Debuging test:
#define NOT_IN_UESED -1
#define DISABLED -1

enum typeP {
    ISOTERM,
    ADIABATIC
};


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

#if BOUNDARIES == PERIODIC
    double boxLimits[2 * DIM];
#endif // PERIODIC
};

#endif // GLOBAL_H
