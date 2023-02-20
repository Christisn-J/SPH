//
// Created by Johannes Martin on 23.09.21.
//

#include "../include/H5.h"

void H5::load(const std::string &file){
    HighFive::File h5file(file, HighFive::File::ReadOnly);

    // read datasets from file
    HighFive::DataSet mass = h5file.getDataSet("/m");
    HighFive::DataSet pos = h5file.getDataSet("/x");
    HighFive::DataSet vel = h5file.getDataSet("/v");
    HighFive::DataSet energy = h5file.getDataSet("/u");

    // read data into containers
    mass.read(m);
    pos.read(x);
    vel.read(v);
    energy.read(u);

    // sanity check
    if (x.size() == v.size() && x.size() == m.size() && x.size() == u.size()){
        nParticles = x.size();
    } else {
        throw std::length_error("Length mismatch between mass, position and/or velocity vectors.");
    }

}

void H5::initialize(Particles &particles){
    std::vector<std::vector<double>>::iterator xit = x.begin();
    std::vector<std::vector<double>>::iterator vit = v.begin();
    std::vector<double>::iterator mit = m.begin();
    std::vector<double>::iterator uit = u.begin();

    int pCounter = 0;

    while (xit != x.end()){
        particles.m[pCounter] = *mit;
        particles.u[pCounter] = *uit;
        particles.x[pCounter] = (*xit)[0];
        particles.vx[pCounter] = (*vit)[0];

        #if DIM >= 2
            particles.y[pCounter] = (*xit)[1];
            particles.vy[pCounter] = (*vit)[1];
        #endif //2D

        #if DIM == 3
                particles.z[pCounter] = (*xit)[2];
                particles.vz[pCounter] = (*vit)[2];
        #endif //3D

        ++xit;
        ++uit;
        ++vit;
        ++mit;
        ++pCounter;
    }
}

void H5::save(const std::string &filename, Particles &particles){
    // TODO:
}
