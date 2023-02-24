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
 * 1: fluid in vacum
 * 2: kevin-helmholtz
 * 3: sedov
**/
#define TESTCASE 0

/** define kind of boundaries should be employed:
 * 0: transparent 
 * 1: periodi
 * 2: reflectiv
**/
#define BOUNDARIES 1

/** define debug level to enable additional output:
 * 0: protforce
**/
#define NNS 1

/** define debug level to enable additional output:
 * 0: iostherm (soundspeed is constant)
**/
#define TYP 0



/// define Courant-Friedrichs-Levy number, should be smaller than 1
#define CFL .4

#endif //PARAMETER_H