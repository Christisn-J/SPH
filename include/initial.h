#include <highfive/H5File.hpp>
#include <vector>

#include "global.h"
#include "Particles.h"

 // containers to be filled from hdf5 file
std::vector<double> m {}, u {};
std::vector<std::vector<double>> r {}, v {};
std::vector<int> matId {};
int numberOfParticles { 0 };

void load(const std::string &file);
void initialize(Particles &particles);