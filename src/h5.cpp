//
// Created by Johannes Martin on 23.09.21.
//

#include "../include/h5.h"

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

    // sanity check
    if (r.size() == v.size() && r.size() == m.size() && r.size() == u.size()){
        nParticles = r.size();
    } else {
        throw std::length_error("Length mismatch between mass, position and/or velocity vectors.");
    }

}

void initialize(Particles &particles){
    std::vector<std::vector<double>>::iterator rit = r.begin();
    std::vector<std::vector<double>>::iterator vit = v.begin();
    std::vector<double>::iterator mit = m.begin();
    std::vector<double>::iterator uit = u.begin();

    int pCounter = 0;

    while (rit != r.end()){
        particles.m[pCounter] = *mit;
        particles.u[pCounter] = *uit;
        particles.x[pCounter] = (*rit)[0];
        particles.vx[pCounter] = (*vit)[0];

        #if DIM >= 2
            particles.y[pCounter] = (*rit)[1];
            particles.vy[pCounter] = (*vit)[1];
        #endif //2D

        #if DIM == 3
                particles.z[pCounter] = (*rit)[2];
                particles.vz[pCounter] = (*vit)[2];
        #endif //3D

        ++rit;
        ++uit;
        ++vit;
        ++mit;
        ++pCounter;
    }
}

// TODO:
void save(const std::string &filename, Particles &particles){
    /*
    // open output file
    HighFive::File h5File {filename, 
                                    HighFive::File::ReadWrite |
                                    HighFive::File::Create |
                                    HighFive::File::Truncate };

    // dimensions for datasets containing vectors
    std::vector<size_t> dataSpaceDims(2);
    dataSpaceDims[0] = std::size_t(N); // number of particles
    dataSpaceDims[1] = DIM;

    // create datasets
    // TODO: Create a h5 object holding all meta data
    HighFive::DataSet rhoDataSet = h5File.createDataSet<double>("/rho", HighFive::DataSpace(NUM));
    HighFive::DataSet mDataSet = h5File.createDataSet<double>("/m", HighFive::DataSpace(NUM));
    HighFive::DataSet uDataSet = h5File.createDataSet<double>("/u", HighFive::DataSpace(NUM));
    HighFive::DataSet posDataSet = h5File.createDataSet<double>("/r", HighFive::DataSpace(dataSpaceDims));
    HighFive::DataSet velDataSet = h5File.createDataSet<double>("/v", HighFive::DataSpace(dataSpaceDims));

    // containers for particle data
    std::vector<double> rhoVec(rho, rho+NUM);
    std::vector<double> mVec(m, m+NUM);
    std::vector<double> uVec(u, u+NUM);

    std::vector<std::vector<double>> posVec(NUM);
    std::vector<std::vector<double>> velVec(NUM);

    std::vector<std::vector<double>> rhoGradVec(NUM);

    // fill containers with data
    std::vector<double> posBuf(DIM);
    std::vector<double> velBuf(DIM);
    std::vector<double> rhoGradBuf(DIM);
    for(int i=0; i<N; ++i){
        // position
        posBuf[0] = x[i];
        #if DIM == 2
            posBuf[1] = y[i];
        #endif // 2D

        #if DIM == 3
                posBuf[2] = z[i];
        #endif // 3D
        posVec[i] = posBuf;

        // velocity
        velBuf[0] = vx[i];
        #if DIM == 2
            velBuf[1] = vy[i];
        #endif // 2D
        #if DIM == 3
                velBuf[2] = vz[i];
        #endif // 3D
        velVec[i] = velBuf;

        // density gradient
        rhoGradBuf[0] = rhoGrad[i][0];
        #if DIM == 2
            rhoGradBuf[1] = rhoGrad[i][1];
        #endif // 2D
        #if DIM == 3
                rhoGradBuf[2] = rhoGrad[i][2];
        #endif // 3D
        rhoGradVec[i] = rhoGradBuf;
    }
    // write data
    rhoDataSet.write(rhoVec);
    mDataSet.write(mVec);
    uDataSet.write(uVec);
    posDataSet.write(posVec);
    velDataSet.write(velVec);
    rhoGradDataSet.write(rhoGradVec);
    */
}
