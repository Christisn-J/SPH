//
// Created by Christian Jetter on 20.02.23
//

#ifndef BOUNDARY_H
#define BOUNDARY_H

#include "parameter.h"
#include "Logger.h"

#define PERIODIC 1
#define REFLECTIV 2

class Domain {

public:
    struct Frame {
        Frame(double *bounds) : minX { bounds[0] }, maxX { bounds[DIM] }
#if DIM >= 2
                              ,minY { bounds[1] }, maxY { bounds[DIM+1] }
#endif // 2D

#if DIM == 3
                              ,minZ { bounds[2] }, maxZ { bounds[DIM+2] }
#endif // 3D
                               {}

        double minX;
        double maxX;
#if DIM >= 2
        double minY;
        double maxY;
#endif // 2D

#if DIM == 3
        double minZ;
        double maxZ;
#endif //3D
    };

    Domain(Frame bounds);
    ~Domain();

    Frame bounds;
};

#endif // BOUNDARY_H