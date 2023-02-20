//
// Created by Christian Jetter on 30.01.23
//

#ifndef H5_H
#define H5_H

#include <highfive/H5File.hpp>
#include <vector>
#include <cmath>

#include "global.h"
#include "Particles.h"

// containers to be filled from hdf5 file
extern std::vector<double> m {}, u {};
extern std::vector<std::vector<double>> r {}, v {};

void load(const std::string &file);
void initialize(Particles &particles);
// TODO:
void save(const std::string &file, Particles &particles);

int nParticles { 0 };

#endif // H5_H