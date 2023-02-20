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

class H5 {
public:
    void load(const std::string &file);
    void initialize(Particles &particles);
    void save(const std::string &file, Particles &particles);

    int getN() const { return nParticles; };

private:
    // containers to be filled from hdf5 file
    std::vector<double> m {}, u {};
    std::vector<std::vector<double>> x {}, v {};
    int nParticles { 0 };

};
#endif // H5_H