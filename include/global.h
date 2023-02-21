//
// Created by Christian Jetter on 23.01.23
//

#ifndef GLOBAL_H
#define GLOBAL_H

#include <parameter.h>
#include <string>

#define NOT_IN_UESED -1

enum typeP {
    ISOTERM,
    ADIABATIC
};

enum typeNN {
    PROTFORCE,
};

struct Configuration {
    std::string initFile;
    std::string outDir;
    double timeStep;
    double timeEnd;
    int storeFrequency;
    int maxInteractions;
    double h; // smoothing lenght
    double gamma; // adiabatic index
};

#endif // GLOBAL_H
