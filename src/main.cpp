//
// test
//

// std libraries c++
#include <iostream>
#include <cmath>

// header only libraries
#include <cxxopts.hpp>

// by Johannes Martin
#include "../include/ConfigParser.h"
#include "../include/Logger.h"
#include "../include/h5.h"

#include "../include/global.h"
#include "../include/Particles.h"
#include "../include/kernel.h"
#include "../include/boundary.h"
#include "../include/lib.h"
#include "../include/tools.h"

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
    ConfigParser::Configuration config;

    config.initFile = read.getVal<std::string>("initFile");
    Logger(INFO) << "    > Initial distribution: " << config.initFile;
    config.outDir = read.getVal<std::string>("outDir");
    Logger(INFO) << "    > Output directory: " << config.outDir;
    config.timeStep = read.getVal<double>("timeStep");
    Logger(INFO) << "    > Time step: " << config.timeStep;
    config.timeEnd = read.getVal<double>("timeEnd");
    Logger(INFO) << "    > End of simulation: " << config.timeEnd;
    config.h = read.getVal<double>("kernelSize");
    Logger(INFO) << "    > Using global kernel size h = " << config.h;
    config.gamma = read.getVal<double>("gamma");
    Logger(INFO) << "    > Adiabatic index for ideal gas EOS gamma = " << config.gamma;

// initialize -------------------------------------------------------------------------------------------------
    Logger(INFO) << "Initializing simulation ...";
    // load inital conditions
    load(config.initFile);
    // instantiation Particles
    Particles sampel = Particles(nParticles);
    // initialize the loaded conditions
    initialize(sampel);

    Logger(INFO) << "    > N = " << sampel.N;
    Logger(INFO) << "... done.";

// simulation ------------------------------------------------------------------------------------------------

    Logger(INFO) << "Starting simulation ...";
    run(config, sampel);
    Logger(INFO) << "... done.";    
        
    return 0;
}