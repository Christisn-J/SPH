//
// Created by Christian Jetter on 30.01.23
//

// std libraries c++
#include <iostream>
#include <cmath>

// header only libraries
#include <cxxopts.hpp>

// by Johannes Martin
#include "../include/ConfigParser.h"
#include "../include/Logger.h"
#include "../include/H5.h"

// by Christian Jetter
#include "../include/global.h"
#include "../include/Particles.h"
#include "../include/Boundary.h"
#include "../include/kernel.h"
#include "../include/tools.h"

#include "../include/lib.h"

// instantiation Logger
structlog LOGCFG = {};

int main(int argc, char *argv[]){
// parser ----------------------------------------------------------------------------------------------------
    cxxopts::Options cmdLineOptions { "mlh",
                                      "Demonstrator for the meshless hydrodynamic simulation methods MFV and MFM" };
    cmdLineOptions.add_options()
            ("c,config", "Path to config file", cxxopts::value<std::string>()->default_value("config.info"))
            ("v,verbose", "More printouts for debugging")
            ("s,silent", "Suppress normal printouts")
            ("h,help", "Show this help");

    auto cmdLineOpts = cmdLineOptions.parse(argc, argv);

    if (cmdLineOpts.count("help")) {
        std::cout << cmdLineOptions.help() << std::endl;
        exit(0);
    }

// logger-----------------------------------------------------------------------------------------------------
    // initialize Logger
    LOGCFG.headers = true;
    LOGCFG.level = cmdLineOpts.count("verbose") ? DEBUG : INFO;

    // avoid overlapping parser arguments 
    if (cmdLineOpts.count("silent")){
        if(cmdLineOpts.count("verbose")){
            throw std::invalid_argument("Command line options -s and -v are incompatible");
        } else {
            LOGCFG.level = WARN;
        }
    }

// config ----------------------------------------------------------------------------------------------------
    Logger(INFO) << "Reading configuration ... ";

    // instantiation ConfigParser
    ConfigParser read { cmdLineOpts["config"].as<std::string>() };
    Configuration config;

    config.initFile = read.getVal<std::string>("initFile");
    Logger(INFO) << "    > Initial distribution: " << config.initFile;
    config.outDir = read.getVal<std::string>("outDir");
    Logger(INFO) << "    > Output directory: " << config.outDir;
    config.timeStep = read.getVal<double>("timeStep");
    Logger(INFO) << "    > Time step: " << config.timeStep;
    config.timeEnd = read.getVal<double>("timeStart");
    Logger(INFO) << "    > Start of simulation: " << config.timeStart;
    config.timeEnd = read.getVal<double>("timeEnd");
    Logger(INFO) << "    > End of simulation: " << config.timeEnd;
    config.storeFrequency = read.getVal<int>("storeFrequency");
    Logger(INFO) << "    > Store data to h5 file every " << config.storeFrequency << " steps";
    config.maxInteractions = read.getVal<int>("maxInteractions");
    Logger(INFO) << "    > Max number of Interactions: " << config.maxInteractions;
    config.h = read.getVal<double>("h");
    Logger(INFO) << "    > Using global kernel size h = " << config.h;
    config.gamma = read.getVal<double>("gamma");
    Logger(INFO) << "    > Adiabatic index for ideal gas EOS gamma = " << config.gamma;

#if BOUNDARIES == PERIODIC
    auto boxLimits = read.getObj("boxLimits");
    config.boxLimits[0] = boxLimits.getVal<double>("lowerX");
    config.boxLimits[DIM] = boxLimits.getVal<double>("upperX");

#if DIM >= 2
    config.boxLimits[1] = boxLimits.getVal<double>("lowerY");
    config.boxLimits[DIM+1] = boxLimits.getVal<double>("upperY");
#endif // 2D

#if DIM == 3
    config.boxLimits[2] = boxLimits.getVal<double>("lowerZ");
    config.boxLimits[DIM+2] = boxLimits.getVal<double>("upperZ");
#endif // 3D

    // TODO: lib/str <double>array2str(<double> arr[]]) 
    std::string ArrayStr = "[";
    for (int i=0; i<2*DIM; i++){
        ArrayStr.append(std::to_string(config.boxLimits[i]));
        if(i<2*DIM-1) ArrayStr.append(", ");
    }
    Logger(INFO) << "    > Periodic boundaries within box: " << ArrayStr << "]";

#endif // PERIODIC

// initialize -------------------------------------------------------------------------------------------------
    H5 distribuition;
    Logger(INFO) << "Initializing simulation ...";
    // load inital conditions
    distribuition.load(config.initFile);
    // instantiation Particles
    Particles sampel = Particles(distribuition.getN(), config);
    // initialize the loaded conditions
    distribuition.initialize(sampel);

    Logger(INFO) << "    > N = " << sampel.N;

#if BOUNDARIES == PERIODIC
    double *domainLimits = config.boxLimits;
#else // PERIODIC
    double domainLimits[DIM*2];
    particles.getDomainLimits(domainLimits);
#endif
    Domain::Frame boundingBox { domainLimits };

    Logger(INFO) << "... done.";

// simulation ------------------------------------------------------------------------------------------------
    Logger(INFO) << "Starting simulation ...";
    algorithm(config, sampel, boundingBox);
    Logger(INFO) << "... done.";    
        
    return 0;
}