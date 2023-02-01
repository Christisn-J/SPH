//
// Created by Johannes Martin on 21.09.22.
//

#ifndef PARAMETER_H
#define PARAMETER_H

/** define debug level to enable additional output:
 * 0: no debug additions
 * 1: additional checks
 * 2: dump NNL and ghosts to files
**/
#define DEBUG_LVL 2

/** define kind of boundaries should be employed:
 * 0: transparent 
 * 1: periodic
 * 2: reflectiv
**/
#define BOUNDARIES 0

/// possible values 1 up to 3 for dimension of the simulations
#define DIM 2

/// define if timestep is adaptive
#define ADAPTIVE_TIMESTEP false

/// define Courant-Friedrichs-Levy number, should be smaller than 1
#define CFL .4

#endif //PARAMETER_H