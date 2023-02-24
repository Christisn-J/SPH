//
// Created by Christian Jetter on 30.01.23
//

#include "../include/tools.h"

void algorithm(Configuration config, Particles particles, Domain::Cell bounds){
    double t = config.timeStart;
    double timeStep;
    int step = 0;

// fixed grid ------------------------------------------------------------------------------------------------
    // instantiation grid
    Domain  domain(bounds);

    Logger(INFO) << "    > Creating grid ... ";
    domain.createGrid(config.h);
    Logger(INFO) << "    > ... got " << domain.numGridCells << " cells";

// start simulation ------------------------------------------------------------------------------------------------
    do {
        Logger(INFO) << "  > TIME: " << t << ", STEP: " << step << " ------------------------------------------------------------------------------------";

#if BOUNDARIES == TRANSPARENT
// changeing grid ------------------------------------------------------------------------------------------------
        Logger(INFO) << "    > Computing domain limits ...";
        double domainLimits[DIM*2];
        particles->getDomainLimits(domainLimits);
        Domain::Cell boundingBox { domainLimits };
        domain.bounds = boundingBox;
        /// TODO: domain.printout();
        Logger(DEBUG) << "      > ... creating grid ...";
        domain.createGrid(config.h);
        Logger(INFO) << "    > ... done.";
#endif // TRANSPARENT

#if NNS == 1
        Logger(INFO) << "    > Assigning particles ...";
        particles.assignParticlesAndCells(domain);
        Logger(INFO) << "    > ... done.";
#endif // NNS GRID

// determine -------------------------------------------------------------------------------------------------
        Logger(INFO) << "    > Nearest neighbor search";
        particles.compNN(domain, config.h);

        Logger(INFO) << "    > Computing density";
        particles.compDensity(config.h);

        Logger(INFO) << "    > Computing pressure";
        particles.compPressure(config.gamma);

// conservation quantities ----------------------------------------------------------------------------------
        Logger(DEBUG) << "      SANITY CHECK > M_tot = " << particles.sumMass();
        Logger(DEBUG) << "      SANITY CHECK > px_tot = " << particles.sumMomentum(1);
#if DIM >= 2
        Logger(DEBUG) << "      SANITY CHECK > py_tot = " << particles.sumMomentum(2);
#endif //2D
#if DIM == 3
        Logger(DEBUG) << "      SANITY CHECK > pz_tot = " << particles.sumMomentum(3);
#endif // 3D
        Logger(DEBUG) << "      SANITY CHECK > E_tot = " << particles.sumEnergy();

// time step -------------------------------------------------------------------------------------------------
#if TYP != ISOTERM
        // change time step
        Logger(INFO) << "    > Selecting global timestep ... ";
        timeStep = particles->compGlobalTimestep(config.gamma, config.kernelSize);
        Logger(INFO) << "    > dt = " << timeStep << " selected.";
#else   
        // constant time step
        timeStep = config.timeStep;
#endif

// gradient---------------------------------------------------------------------------------------------------
        Logger(INFO) << "    > Computing gradients";
        particles.gradient(particles.rho, particles.rhoGrad);
        particles.gradient(particles.vx, particles.vxGrad);
#if DIM >= 2
        particles.gradient(particles.vy, particles.vyGrad);
#endif // 2D
#if DIM == 3
        particles->gradient(particles->vz, particles->vzGrad);
#endif // 3D
        particles.gradient(particles.P, particles.PGrad);

// save dataset ----------------------------------------------------------------------------------------------
      if (step % config.storeFrequency == 0) {
            std::stringstream stepss;
            Logger(INFO) << "   > Dump particle distribution";
            stepss << std::setw(6) << std::setfill('0') << step;
            Logger(INFO) << "      > save particles to file";
            particles.save(config.outDir + "/" + stepss.str() + std::string(".h5"), t);
        }

// break condition -------------------------------------------------------------------------------------------
        if (t>=config.timeEnd){
            Logger(INFO) << "    > t = " << t << " -> FINISHED!";
            break;
        }

// calculation ----------------------------------------------------------------------------------------------- 
        Logger(INFO) << "    > Calculate acceleration";
        particles.accelerate(config.h);
        particles.damping(config.h);
#if BOUNDARIES != TRANSPARENT
// boundary--- ----------------------------------------------------------------------------------------------- 

#endif

// update ----------------------------------------------------------------------------------------------------
        Logger(INFO) << "    > Updating state";     
        particles.integrate(timeStep, domain);

        t += timeStep;
        ++step;

// break condition -------------------------------------------------------------------------------------------
    } while(t<config.timeEnd+timeStep);
}