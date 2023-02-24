//
// Created by Christian Jetter on 02.02.23
//

#include "../include/Particles.h"

// constructor ----------------------------------------------------------------------------------------------- 
Particles::Particles(int nParticles, Configuration config) : N { nParticles }{
    // prevent stack overflow
    if (config.maxInteractions < N){
        MAX_INTERACTIONS = config.maxInteractions;
    }else{
        MAX_INTERACTIONS = N;
    }
    
    // allocate memory
    m = new double[N];          // mass
    u = new double[N];          // energy

    rho = new double[N];        // density
    rhoGrad = new double[N][DIM];

    P = new double[N];          // pressure
    PGrad = new double[N][DIM];

    x = new double[N];          // possition
    vx = new double[N];         // velossity
    vxDelta = new double[N];
    vxGrad = new double[N][DIM];
    
#if DIM >= 2
    y = new double[N];
    vy = new double[N];
    vyDelta = new double[N];
    vyGrad = new double[N][DIM];
#endif // 2D
#if DIM == 3
    z = new double[N];
    vz = new double[N];
    vzDelta = new double[N];
    vzGrad = new double[N][DIM];
#endif // 3D

    nnl = new int[N*MAX_INTERACTIONS];
    noi = new int[N];
    cell = new int[N];

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
    delete[] vxDelta;
    delete[] vxGrad;

#if DIM >= 2
    delete[] y;
    delete[] vy;
    delete[] vyDelta;
    delete[] vyGrad;
#endif // 2D
#if DIM == 3
        delete[] z;
        delete[] vz;
        delete[] vzDelta;
        delete[] vzGrad;
#endif // 3D

    delete[] nnl;
    delete[] noi;
    delete[] cell;

#if BOUNDARIES != TRANSPARENT
#endif
}

// functions -------------------------------------------------------------------------------------------------
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
        // initialiase 
        vxDelta[i] = 0;
#if DIM >= 2
        vyDelta[i] = 0;
#endif // 2D
#if DIM == 3
        vzDelta[i] = 0;
#endif // 3D

#if TESTCASE == DEBUG
        Logger(DEBUG) << i << "\ta " << vxDelta[i] << " / " << vyDelta[i] ;      
#endif // TESTCASE
        for(int n=0; n<noi[i]; ++n){
            int j = nnl[getNNLidx(i,n)];
            double r = distance(i,j);

            if (r == 0 || i == j){
                // no total share because of r := |r_i-r_n| multiplier
                continue;
            }

            // ... dW(r)/r * (x_i-x_n) with r := r_i-r_n
            vxDelta[i] += -m[j]*(P[i]/pow(rho[i],2)+P[j]/pow(rho[j],2))*dW(r,h)/r* (x[i]-x[j]);

#if DIM >= 2
            vyDelta[i] += -m[j]*(P[i]/pow(rho[i],2)+P[j]/pow(rho[j],2))*dW(r,h)/r* (y[i]-y[j]);
#endif // 2D
#if DIM == 3
            vzDelta[i] += -m[j]*(P[i]/pow(rho[i],2)+P[j]/pow(rho[j],2))*dW(r,h)/r* (z[i]-z[j]);
#endif // 3D
        }

#if TESTCASE == DEBUG
        Logger(DEBUG) << i << "\ta " << vxDelta[i] << " / " << vyDelta[i] ;      
#endif // TESTCASE
    }
}

void Particles::damping(const double &h){
    for(int i=0; i<N; ++i){

#if TESTCASE == DEBUG
        vxDelta[i] += -vxDelta[i];
#if DIM >= 2
        vyDelta[i] += -vyDelta[i];
#endif // 2D
#if DIM == 3
        vzDelta[i] += -vzDelta[i];
#endif // 3D

#endif // TESTCASE
    }
}

void Particles::integrate(const double &dt, const Domain &domain){
    // Eulerstep
    // update velossity
    for(int i=0; i<N; ++i ){
#if TESTCASE == DEBUG
        Logger(DEBUG) << i << "\tv_n\t" << vx[i] << " / " << vy[i];
#endif // TESTCASE

        vx[i] = vx[i] + dt * vxDelta[i]; 
#if DIM >= 2
        vy[i] = vy[i] + dt * vyDelta[i];
#endif // 2D
#if DIM == 3
        vz[i] = vz[i] + dt * vzDelta[i];
#endif // 3D

#if TESTCASE == DEBUG
        Logger(DEBUG) << "\tv_n+1\t" << vx[i] << " / " << vy[i];
#endif // TESTCASE

    }

    // update possition
    for(int i=0; i<N; ++i ){
#if TESTCASE == DEBUG
        Logger(DEBUG) << i << "\tpos_n\t" << x[i] << " / " << y[i] ;
#endif // TESTCASE

        x[i] = x[i]+ dt*vx[i];
#if DIM >= 2
        y[i] = y[i]+ dt*vy[i];
#endif // 2D
#if DIM == 3
        z[i] = z[i]+ dt*vz[i];
#endif // 3D

#if TESTCASE == DEBUG
        Logger(DEBUG) << "\tpos_n+1\t" << x[i] << " / " << y[i] ;
#endif // TESTCASE

#if BOUNDARIES == PERIODIC
        if (x[i] < domain.bounds.minX) {
            x[i] = domain.bounds.maxX - (domain.bounds.minX - x[i]);
        } else if (domain.bounds.maxX <= x[i]) {
            x[i] = domain.bounds.minX + (x[i] - domain.bounds.maxX);
        }
#if DIM >= 2
        if (y[i] < domain.bounds.minY) {
            y[i] = domain.bounds.maxY - (domain.bounds.minY - y[i]);
        } else if (domain.bounds.maxY <= y[i]) {
            y[i] = domain.bounds.minY + (y[i] - domain.bounds.maxY);
        }
#endif // 2D
#if DIM ==3
        if (z[i] < domain.bounds.minZ) {
            z[i] = domain.bounds.maxZ - (domain.bounds.minZ - z[i]);
        } else if (domain.bounds.maxZ <= z[i]) {
            z[i] = domain.bounds.minZ + (z[i] - domain.bounds.maxZ);
        }
#endif // 3D
#endif // PERIODIC

    }
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
#if DIM >= 2
        posBuf[1] = y[i];
#endif // 2D
#if DIM == 3
        posBuf[2] = z[i];
#endif // 3D
        posVec[i] = posBuf;

        // velocity
        velBuf[0] = vx[i];
#if DIM >= 2
        velBuf[1] = vy[i];
#endif // 2D
#if DIM == 3
        velBuf[2] = vz[i];
#endif // 3D
        velVec[i] = velBuf;

        // density gradient
        rhoGradBuf[0] = rhoGrad[i][0];
#if DIM >= 2
        rhoGradBuf[1] = rhoGrad[i][1];
#endif // 2D
#if DIM == 3
        rhoGradBuf[2] = rhoGrad[i][2];
#endif // 3D
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

// NNS functions -------------------------------------------------------------------------------------------------------
#if NNS == GRID
void Particles::assignParticlesAndCells(Domain &domain){

    // reset particles assigned to grid cells
    for(int iGrid=0; iGrid<domain.numGridCells; ++iGrid){
        domain.grid[iGrid].prtcls = std::vector<int>();
    }

    for(int i=0; i<N; ++i){

        // This rarely happens when x or y is really close to the domain bounds
        int floorX = floor((x[i]-domain.bounds.minX)/domain.cellSizeX);
        if (floorX == domain.cellsX){
            floorX -= 1;
        }
#if DIM >= 2
        int floorY = floor((y[i]-domain.bounds.minY)/domain.cellSizeY);
        if (floorY == domain.cellsY){
            floorY -= 1;
        }
#endif // 2D
#if DIM == 3
        int floorZ = floor((z[i]-domain.bounds.minZ)/domain.cellSizeZ);
        if (floorZ == domain.cellsZ){
            floorZ -= 1;
        }
#endif // 3D


        int iGrid = floorX 
#if DIM >= 2
        + floorY * domain.cellsX
#endif // 2D
#if DIM == 3
        + floorZ * domain.cellsX * domain.cellsY
#endif // 3D
        ;
        
        domain.grid[iGrid].prtcls.push_back(i);
        cell[i] = iGrid; // assign cells to particles
    }
}
#endif // NNS GRID

void Particles::compNN(Domain &domain, const double &h){
    // loop over particles
    for(int i=0; i<N; ++i){
        int interact = 0;

#if NNS == PROTFORCE
        // seachr nearest neighbor
        for(int n=0; n<N; ++n){
            if(i == n){
                // not interact with each self
                continue;
            }

            // cheack distance beetween particles
            if(distance(i,n) < h){
                if (interact >= N) {
                    Logger(ERROR) << "Out of bounds: Number of particles";
                    exit(1);
                }
                if (interact >= MAX_INTERACTIONS) {
                    Logger(ERROR) << "To many interactions";
                    exit(1);
                }

                nnl[i * MAX_INTERACTIONS+interact] = n;
                ++ interact;
                
            }
        }
        noi[i] = interact;
#if TESTCASE == DEBUG      
        Logger(DEBUG) << i  << "\tnoi "<< noi[i];
#endif // TESTCASE
#endif // NNS PROTFORCE

#if NNS == GRID
        int numSearchCells = pow(3, DIM);
        // search for nearest neighbors in the particle cell and neighbor cells
        int cells[numSearchCells];
        // neighboring cells (including particles cell)
        domain.getNeighborCells(cell[i], cells);
        // do nearest neighbor search
        for (int iNeighbor=0; iNeighbor<numSearchCells; ++iNeighbor){
            // loop over particle indices in all
            if (cells[iNeighbor] < 0){
                /// TODO: handle ghost cells in external function
            } else {
                for(auto const &iPrtcl : domain.grid[cells[iNeighbor]].prtcls){
                    if (iPrtcl != i){
                        if (distance(i,iPrtcl)  < h) {
                            if (interact >= MAX_INTERACTIONS) {
                                Logger(ERROR) << "MAX_NUM_INTERACTIONS exceeded for particle "
                                              << i << " - Aborting.";
                                exit(1);
                            }
                            nnl[i * MAX_INTERACTIONS+interact] = iPrtcl;
                            ++ interact;
                        }
                    }
                }
            }
        }
        if (interact == 0){
            Logger(WARN) << "No neighbors for particle " << i << ". Caution.";
        }
        noi[i] = interact;
#endif // NNS GRID

    }
}

// boundary functions-------------------------------------------------------------------------------------------------------
#if BOUNDARIES == TRANSPARENT
void Particles::getDomainLimits(double *domainLimits){

    double minX { std::numeric_limits<double>::max() };
    double maxX { std::numeric_limits<double>::min() };
#if DIM >= 2
    double minY { std::numeric_limits<double>::max() };
    double maxY { std::numeric_limits<double>::min() };
#endif // 2D
#if DIM == 3
    double minZ { std::numeric_limits<double>::max() };
    double maxZ { std::numeric_limits<double>::min() };
#endif // 3D

    for(int i=0; i<N; ++i){
        if (x[i] < minX){
            minX = x[i];
        } else if (x[i] > maxX){
            maxX = x[i];
        }
#if DIM >= 2
        if (y[i] < minY){
            minY = y[i];
        } else if (y[i] > maxY){
            maxY = y[i];
        }
#endif // 2D
#if DIM == 3
        if (z[i] < minZ){
            minZ = z[i];
        } else if (z[i] > maxZ){
            maxZ = z[i];
        }
#endif // 3D
    }

    domainLimits[0] = minX;
    domainLimits[DIM] = maxX;
#if DIM >= 2
    domainLimits[1] = minY;
    domainLimits[DIM+1] = maxY;
#endif // 2D
#if DIM == 3
    domainLimits[2] = minZ;
    domainLimits[DIM+2] = maxZ;
#endif // 3D
}
#endif // TRANSPARENT

// compelation functions -------------------------------------------------------------------------------------------------------
void Particles::compDensity(const double &h){
    for(int i=0; i<N; ++i){
        // initialiase
        rho[i] = 0; 

        for(int n=0; n<noi[i]; ++n){
            int j = nnl[getNNLidx(i,n)];

            rho[i] += m[j]*W(distance(i,j),h);
        }

#if TESTCASE == DEBUG
            Logger(DEBUG) << i << "\trho "<< rho[i];      
#endif // TESTCASE
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
#else // ISOZTERM
        P[i] = (gamma-1.)*rho[i]*u[i];
#endif // NOT ISOTERM

        if(P[i] <= 0.){
            Logger(WARN) << "Zero or negative pressure @" << i;
        }
    }
}

// sanity check functions ------------------------------------------------------------------------------------
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

// helper functions ------------------------------------------------------------------------------------------
int Particles::getNNLidx(int i, int n){
    return i*MAX_INTERACTIONS+n;
}

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
#endif// 3D