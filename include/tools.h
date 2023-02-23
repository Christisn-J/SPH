//
// Created by Christian Jetter on 23.01.23
//

#ifndef TOOLS_H
#define TOOLS_H

#include "parameter.h"
#include "Logger.h"
#include "Particles.h"
#include "ConfigParser.h"
#include "Boundary.h"
#include "lib.h"

void algorithm(Configuration config, Particles particles, Domain::Frame bounds);

#endif // TOOLS_H