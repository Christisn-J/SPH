//
// Created by Christian Jetter on 02.02.23
//

#include "../include/Particles.h"

// constructor ----------------------------------------------------------------------------------------------- 
Particles::Particles(int nParticles, Configuration config) : N { nParticles }, MAX_INTERACTIONS{ config.maxInteractions }{
    // allocate memory
    m = new double[N];
    u = new double[N];

    rho = new double[N];
    rhoGrad = new double[N][DIM];

    P = new double[N];
    PGrad = new double[N][DIM];

    x = new double[N];
    vx = new double[N];
    vxGrad = new double[N][DIM];
    
#if DIM >= 2
    y = new double[N];
    vy = new double[N];
    vyGrad = new double[N][DIM];
#endif // 2D

#if DIM == 3
    z = new double[N];
    vz = new double[N];
    vzGrad = new double[N][DIM];
#endif // 3D

    nnl = new int[N*MAX_INTERACTIONS];
    noi = new int[N];
}

// deconstructor --------------------------------------------------------------------------------------------- 
Particles::~Particles() {
    delete[] m;
    delete[] u;

    delete[] rho;
    delete[] rhoGrad;

    delete[] P;
    delete[] PGrad;

    delete[] x;
    delete[] vx;
    delete[] vxGrad;

#if DIM >= 2
    delete[] y;
    delete[] vy;
    delete[] vyGrad;
#endif // 2D

#if DIM == 3
        delete[] z;
        delete[] vz;
        delete[] vzGrad;
#endif // 3D

#if BOUNDARIES == PERIODIC
#endif
}

// functions -------------------------------------------------------------------------------------------------
void Particles::compNN(const double &h){
    // loop over particles
    for(int i=0; i<N; ++i){
#if NN_SEARCH == PROTFORCE
        // seachr nearest neighbor
        int interact = 0;
        
        for(int n=0; n<N; ++n){
            // cheack distance beetween particles
            if(distance(i,n) < h){
                if (interact >= N) {
                    Logger(ERROR) << "Out of bounds: Number of particles";
                    exit(1);
                    }
                
                nnl[i * N+interact] = n;
                ++ interact;
                
            }
        }
        noi[i] = interact;
        
        Logger(DEBUG) << i  << " noi "<< noi[i];
#endif // PROTFORCE
    }

}

void Particles::compDensity(const double &h){
    for(int i=0; i<N; ++i){
        rho[i] = m[i]*kernel(distance(i,nnl[i]),h);
#if DEBUG_LVL == 2
        Logger(DEBUG) << i << " rho "<< rho[i]; 
#endif
        if(rho[i] <= 0.){
            Logger(WARN) << "Zero or negative density @" << i;
        }
    }
}

void Particles::compPressure(const double &gamma){
    for (int i=0; i<N; ++i){

#if TYP == ISOTERM
        double c_s = 1.0;
        P[i] = pow(c_s,2)* rho[i];
#endif

#if TYP != ISOTERM
        P[i] = (gamma-1.)*rho[i]*u[i];
#endif

        if(P[i] <= 0.){
            Logger(WARN) << "Zero or negative pressure @" << i;
        }
    }
}

// Sanity check functions ------------------------------------------------------------------------------------
double Particles::sumMass(){
    double M = 0.;
    for (int i=0; i<N; ++i){
        if (std::isnan(m[i])){
            Logger(WARN) << "!! m[" << i <<"] " << "is nan. !!";
        }
        M += m[i];
    }
    return M;
}

double Particles::sumMomentum(const int flag){
    switch (flag){
        case 1:
            return sumMomentumX();
            break;

#if DIM >= 2
        case 2:
            return sumMomentumY();
            break;
#endif // 2D

#if DIM == 3
        case 3:
            return sumMomentumZ();
            break;
#endif // 3D
        default:
            Logger(ERROR) << "No such high Dimension is supported.";
            exit(1);
    }
}

double Particles::sumEnergy(){
    double E = 0.;
    for (int i=0; i<N; ++i){
#if DIM == 1
        E += m[i]*(u[i] + .5*(vx[i]*vx[i]));
#endif // 1D

#if DIM == 2
        E += m[i]*(u[i] + .5*(vx[i]*vx[i]+vy[i]*vy[i]));
#endif // 2D

#if DIM == 3
        E += m[i]*(u[i] + .5*(vx[i]*vx[i]+vy[i]*vy[i]+ vz[i]*vz[i]));
#endif // 3D
    }
    return E;
}

// -----------------------------------------------------------------------------------------------------------
double Particles::distance(const int i, const int n){
    double sum = 0; 
#if DIM >= 1
    sum += pow(x[i]-x[n],2);
#endif // 1D

#if DIM >= 2
    sum += pow(y[i]-y[n],2);
#endif // 2D

#if DIM == 3
    sum += pow(z[i]-z[n],2);
#endif // 3D

    return sqrt(sum);
}

void Particles::gradient(double *f, double (*grad)[DIM]){
    for (int i=0; i<N; ++i) {
        // set values to zerros
        for (int dim = 0; dim < DIM; ++dim) {
            grad[i][dim] = 0;
        }

        // calculate gardiant between NN
        for (int j = 0; j < noi[i]; ++j) {
            for (int dim = 0; dim < DIM; ++dim) {
                grad[i][dim] += (f[nnl[j + i]] - f[i]);
            }
        }
    }
}

void Particles::accelerate(const double &h){
    for(int i=0; i<N; ++i){
        for(int n=0; n<nnl[i]; n++){
            vx[i] += -m[i]*(P[i]/pow(rho[i],2)-P[n]/pow(rho[n],2))*dKernel(distance(i,nnl[i]),h);
        
#if DIM >= 2
            vy[i] += -m[i]*(P[i]/pow(rho[i],2)-P[n]/pow(rho[n],2))*dKernel(distance(i,nnl[i]),h);
#endif // 2D

#if DIM == 3
            vz[i] += -m[i]*(P[i]/pow(rho[i],2)-P[n]/pow(rho[n],2))*dKernel(distance(i,nnl[i]),h);
#endif // 3D
        }
    }
}

void Particles::damping(const double &h){
    for(int i=0; i<N; ++i){
        x[i] = x[i];
        vx[i] = vx[i];

#if DIM >= 2
        y[i] = y[i];
        vy[i] = vy[i];
#endif // 2D

#if DIM == 3
        z[i] = z[i];
        vz[i] = vz[i];
#endif // 3D
    }
}

void Particles::integrate(const double &t){
    // TODO
}

void Particles::save(std::string filename, double simTime){
    // open output file
    HighFive::File h5File { filename, HighFive::File::ReadWrite |
                                      HighFive::File::Create |
                                      HighFive::File::Truncate };

    // dimensions for datasets containing vectors
    std::vector<size_t> dataSpaceDims(2);
    dataSpaceDims[0] = std::size_t(N); // number of particles
    dataSpaceDims[1] = DIM;

    // create datasets
    // TODO: Create a h5 object holding all meta data
    HighFive::DataSet timeDataSet = h5File.createDataSet<double>("/time", HighFive::DataSpace(1));
    HighFive::DataSet rhoDataSet = h5File.createDataSet<double>("/rho", HighFive::DataSpace(N));
    HighFive::DataSet mDataSet = h5File.createDataSet<double>("/m", HighFive::DataSpace(N));
    HighFive::DataSet uDataSet = h5File.createDataSet<double>("/u", HighFive::DataSpace(N));
    HighFive::DataSet posDataSet = h5File.createDataSet<double>("/x", HighFive::DataSpace(dataSpaceDims));
    HighFive::DataSet velDataSet = h5File.createDataSet<double>("/v", HighFive::DataSpace(dataSpaceDims));
    HighFive::DataSet rhoGradDataSet = h5File.createDataSet<double>("/rhoGrad", HighFive::DataSpace(dataSpaceDims));
    HighFive::DataSet PDataSet = h5File.createDataSet<double>("/P", HighFive::DataSpace(N));
    HighFive::DataSet noiDataSet = h5File.createDataSet<int>("/noi", HighFive::DataSpace(N));


    // containers for particle data
    std::vector<double> timeVec({ simTime });
    std::vector<double> rhoVec(rho, rho+N);
    std::vector<double> mVec(m, m+N);
    std::vector<double> uVec(u, u+N);
    std::vector<double> PVec(P, P+N);
    std::vector<int> noiVec(noi, noi+N);

    std::vector<std::vector<double>> posVec(N);
    std::vector<std::vector<double>> velVec(N);

    std::vector<std::vector<double>> rhoGradVec(N);

    // fill containers with data
    std::vector<double> posBuf(DIM);
    std::vector<double> velBuf(DIM);
    std::vector<double> rhoGradBuf(DIM);
    for(int i=0; i<N; ++i){
        //Logger(DEBUG) << "      > Dumping particle @"  << i;
        // position
        posBuf[0] = x[i];
        posBuf[1] = y[i];
#if DIM == 3
        posBuf[2] = z[i];
#endif
        posVec[i] = posBuf;

        // velocity
        velBuf[0] = vx[i];
        velBuf[1] = vy[i];
#if DIM == 3
        velBuf[2] = vz[i];
#endif
        velVec[i] = velBuf;

        // density gradient
        rhoGradBuf[0] = rhoGrad[i][0];
        rhoGradBuf[1] = rhoGrad[i][1];
#if DIM == 3
        rhoGradBuf[2] = rhoGrad[i][2];
#endif
        rhoGradVec[i] = rhoGradBuf;
    }
    // write data
    timeDataSet.write(timeVec); // dummy vec containing one element
    rhoDataSet.write(rhoVec);
    mDataSet.write(mVec);
    uDataSet.write(uVec);
    PDataSet.write(PVec);
    noiDataSet.write(noi);
    posDataSet.write(posVec);
    velDataSet.write(velVec);
    rhoGradDataSet.write(rhoGradVec);
}

// helper functions ------------------------------------------------------------------------------------------
double Particles::sumMomentumX(){
    double momX = 0.;
    for (int i=0; i<N; ++i){
        momX += m[i]*vx[i];
    }
    return momX;
}

#if DIM >= 2
double Particles::sumMomentumY(){
    double momY = 0.;
    for (int i=0; i<N; ++i){
        momY += m[i]*vy[i];
    }
    return momY;
}
#endif // 2D

#if DIM == 3
double Particles::sumMomentumZ(){
    double momZ = 0.;
    for (int i=0; i<N; ++i){
        momZ += m[i]*vz[i];
    }
    return momZ;
}
#endif// 2D