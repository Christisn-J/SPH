//
// Created by Johannes Martin on 21.09.22.
//

#ifndef DOMAIN_H
#define DOMAIN_H

#include <cmath>
#include <vector>

#include "global.h"
#include "Logger.h"

class Domain {

public:
    struct Cell {
        Cell(double *bounds) : minX { bounds[0] }, maxX { bounds[DIM] }
#if DIM >= 2
                                ,minY { bounds[1] }, maxY { bounds[DIM+1] }
#endif // 2D
#if DIM == 3
                                ,minZ { bounds[2] }, maxZ { bounds[DIM+2] }
#endif // 3D
                               {}

        Cell() : minX { 0. }, maxX { 0. }
#if DIM >= 2        
                ,minY { 0. }, maxY { 0. }
#endif // 2D
#if DIM == 3
                ,minZ { 0. }, maxZ { 0. }
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
#endif // 3D
        std::vector<int> prtcls {};
    };

    Domain(Cell bounds);
    ~Domain();

    int numGridCells { 0 };
    std::vector<Cell> grid;
    void createGrid(const double &kernelSize);
    void getNeighborCells(const int &iCell, int *neighborCells);

    Cell bounds; // global cell

// number of grid cells in each dimension
    int cellsX { 0 };
#if DIM >= 2
    int cellsY { 0 };
#endif // 2D
#if DIM == 3
    int cellsZ { 0 };
#endif // 3D

   // number of grid cells in each dimension
    double cellSizeX { 0 };
#if DIM >= 2
    double cellSizeY { 0 };
#endif // 2D
#if DIM == 3
    double cellSizeZ { 0 };
#endif // 3D

private:
    int (*dimIndex)[DIM] { nullptr };    
};


#endif // DOMAIN_H
