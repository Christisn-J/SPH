#include "initial.h"

void load(const std::string &file){
    HighFive::File h5file(file, HighFive::File::ReadOnly);


    // read datasets from file
    HighFive::DataSet mass = h5file.getDataSet("/m");
    HighFive::DataSet pos = h5file.getDataSet("/r");
    HighFive::DataSet vel = h5file.getDataSet("/v");
    HighFive::DataSet materialId = h5file.getDataSet("/materialId");
    HighFive::DataSet energy = h5file.getDataSet("/u");

    // read data into containers
    mass.read(m);
    pos.read(r);
    vel.read(v);
    energy.read(u);
    materialId.read(matId);

}

void initialize(Particles &particles){
    std::vector<std::vector<double>>::iterator rit = r.begin();
    std::vector<std::vector<double>>::iterator vit = v.begin();
    std::vector<double>::iterator mit = m.begin();
    std::vector<double>::iterator uit = u.begin();
    std::vector<int>::iterator matIdIt = matId.begin();

    int pCounter = 0;

    while (rit != x.end()){
        particles.m[pCounter] = *mit;
        particles.u[pCounter] = *uit;
        particles.matId[pCounter] = *matIdIt;
        particles.x[pCounter] = (*rit)[0];
        particles.vx[pCounter] = (*vit)[0];

        #if DIM == 2
        particles.y[pCounter] = (*rit)[1];
        particles.vy[pCounter] = (*vit)[1];
        #endif //2D

        #if DIM == 3
                particles.z[pCounter] = (*rit)[2];
                particles.vz[pCounter] = (*vit)[2];
        #endif //3D

        ++rit;
        ++matIdIt;
        ++uit;
        ++vit;
        ++mit;
        ++pCounter;
    }
}