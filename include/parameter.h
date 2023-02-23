//
// Created by Johannes Martin on 21.09.22.
//

#ifndef PARAMETER_H
#define PARAMETER_H

/// possible values 1 up to 3 for dimension of the simulations
#define DIM 2

/** define debug level to enable additional output:
 * 0: no debug additions
 * 1: additional checks
 * 2: dump NNL and ghosts to files
**/
#define DEBUG_LVL 2

/** define debug level to enable additional output:
 * 0: fluid in vacum
 * 1: fluid with viscosity (toystar)
 * 2: kevin-helmholtz
 * 3: sedov
**/
#define TESTCASE -1

/** define kind of boundaries should be employed:
 * 0: transparent 
 * 1: periodic
 * 2: reflectiv
**/
#define BOUNDARIES 1

/** define debug level to enable additional output:
 * 0: protforce
**/
#define NN_SEARCH 0

/** define debug level to enable additional output:
 * 0: iostherm (soundspeed is constant)
**/
#define TYP 0



/// define Courant-Friedrichs-Levy number, should be smaller than 1
#define CFL .4

#endif //PARAMETER_H