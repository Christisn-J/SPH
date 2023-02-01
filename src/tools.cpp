//
// Created by Christian Jetter on 30.01.23
//

#include "../include/tools.h"

void run(ConfigParser::Configuration config, Particles particles){
    double t = 0;
    double timeStep;
    int step = 0;

    do {
        Logger(INFO) << "  > TIME: " << t << ", STEP: " << step;

// determine -------------------------------------------------------------------------------------------------
        Logger(INFO) << "    > Nearest neighbor search";
        particles.compNN(config.h);

        Logger(INFO) << "    > Computing density";
        particles.compDensity(config.h);

        Logger(INFO) << "    > Computing pressure";
        particles.compPressure(config.gamma);

        Logger(DEBUG) << "      SANITY CHECK > V_tot = " << particles.sumVolume();
        Logger(DEBUG) << "      SANITY CHECK > M_tot = " << particles.sumMass();
        Logger(DEBUG) << "      SANITY CHECK > E_tot = " << particles.sumEnergy();
        Logger(DEBUG) << "      SANITY CHECK > px_tot = " << particles.sumMomentumX();
#if DIM >= 2
        Logger(DEBUG) << "      SANITY CHECK > py_tot = " << particles.sumMomentumY();
#endif //2D

#if DIM == 3
        Logger(DEBUG) << "      SANITY CHECK > pz_tot = " << particles.sumMomentumZ();
#endif // 3D

// time step -------------------------------------------------------------------------------------------------
#if ADAPTIVE_TIMESTEP
        Logger(INFO) << "    > Selecting global timestep ... ";
        timeStep = particles->compGlobalTimestep(config.gamma, config.kernelSize);
        Logger(INFO) << "    > dt = " << timeStep << " selected.";
#else
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
      // TODO:
      if (step % config.dumpInterval == 0) {
            std::stringstream stepss;
            Logger(INFO) << "   > Dump particle distribution";
            stepss << std::setw(6) << std::setfill('0') << step;
            Logger(INFO) << "      > save particles to file";
            save(config.outDir + "/" + stepss.str() + std::string(".h5"), t);
        }

// break condition -------------------------------------------------------------------------------------------
        if (t>=config.timeEnd){
            Logger(INFO) << "    > t = " << t << " -> FINISHED!";
            break;
        }

// calculation ----------------------------------------------------------------------------------------------- 
        Logger(INFO) << "    > Calculate acceleration";
        particles.accelerate(particles.rho, particles.P)
        particles.accelerate(particles.x)
        particles.damping(particles.vx)

#if DIM >= 2
        particles.accelerate(particles.y)
        particles.damping(particles.vy)
#endif // 2D

#if DIM == 3
        particles.accelerate(particles.z)
        particles.damping(particles.vz)
#endif // 3D
    
    


// update ----------------------------------------------------------------------------------------------------
        Logger(INFO) << "    > Updating state";
        particles.integrator(timeStep);

        t += timeStep;
        ++step;

    } while(t<config.timeEnd+timeStep);
}