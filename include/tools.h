//
// Created by Christian Jetter on 23.01.23
//

#ifndef TOOLS_H
#define TOOLS_H

#include "global.h"
#include "ConfigParser.h"
#include "Particles.h"
#include "Logger.h"
#include "Domain.h"

void algorithm(Configuration config, Particles particles, Domain::Cell bounds);

#endif // TOOLS_H