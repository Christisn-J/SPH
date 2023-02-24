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
<<<<<<< HEAD
 * 1: fluid in vacum
=======
 * 0: fluid in vacum
 * 1: TODO: fluid with viscosity (toystar)    
>>>>>>> 8fb93e9614a3af6b49b8ffa2fd130fb67c347e2c
 * 2: kevin-helmholtz
 * 3: TODO: sedov
**/
#define TESTCASE 0

/** define kind of boundaries should be employed:
 * 0: transparent 
<<<<<<< HEAD
 * 1: periodi
 * 2: reflectiv
=======
 * 1: periodic
 * 2: TODO: reflectiv
>>>>>>> 8fb93e9614a3af6b49b8ffa2fd130fb67c347e2c
**/
#define BOUNDARIES 1

/** define debug level to enable additional output:
 * 0: protforce
 * 1: TODO: ---
**/
#define NNS 1

/** define debug level to enable additional output:
 * 0: iostherm (soundspeed is constant)
 * 1: TODO: ---
**/
#define TYP 0

/// define Courant-Friedrichs-Levy number, should be smaller than 1
#define CFL .4

#endif //PARAMETER_H