//
// Created by Johannes Martin on 21.09.22.
//

#include "../include/Domain.h"

Domain::Domain(Cell bounds): bounds { bounds }{}

Domain::~Domain(){
    delete[] dimIndex;
}

void Domain::createGrid(const double &kernelSize){
    cellsX = floor((bounds.maxX - bounds.minX)/kernelSize);
#if DIM >= 2
    cellsY = floor((bounds.maxY - bounds.minY)/kernelSize);
#endif // 2D
#if DIM == 3
    cellsZ = floor((bounds.maxZ - bounds.minZ)/kernelSize);
#endif

    numGridCells = cellsX;
#if DIM >= 2
    numGridCells *= cellsY;
#endif // 2D
#if DIM == 3
    numGridCells *= ellsZ;
#endif

    cellSizeX = (bounds.maxX - bounds.minX)/(double)cellsX;
    Logger(DEBUG) << "      > cellSizeX = " << cellSizeX << ", cellsX = " << cellsX;
#if DIM >= 2
    cellSizeY = (bounds.maxY - bounds.minY)/(double)cellsY;
    Logger(DEBUG) << "      > cellSizeY = " << cellSizeY<< ", cellsY = " << cellsY;
#endif // 2D
#if DIM == 3
    cellSizeZ = (bounds.maxZ - bounds.minZ)/(double)cellsZ;
    Logger(DEBUG) << "      > cellSizeZ = " << cellSizeZ << ", cellsZ = " << cellsZ;
#endif // 3D

    grid = std::vector<Cell>(numGridCells);
    dimIndex = new int[numGridCells][DIM];


    for(int iX=0; iX<cellsX; ++iX){
#if DIM == 1 
        double cellBounds[] = { iX*cellSizeX+bounds.minX, (iX+1)*cellSizeX+bounds.minX};
        grid[iX] = Cell(cellBounds);
        dimIndex[iX][0] = iX;  
#endif // 1D
#if DIM == 2
        for(int iY=0; iY<cellsY; ++iY){
            double cellBounds[] = { iX*cellSizeX+bounds.minX, iY*cellSizeY+bounds.minY,
                                    (iX+1)*cellSizeX+bounds.minX, (iY+1)*cellSizeY+bounds.minY };
            grid[iX+iY*cellsX] = Cell(cellBounds);
            dimIndex[iX+iY*cellsX][0] = iX;
            dimIndex[iX+iY*cellsX][1] = iY;
#endif // 2D
#if DIM == 3
            for(int iZ=0; iZ<cellsZ; ++iZ){
                double cellBounds[] = { iX*cellSizeX+bounds.minX, iY*cellSizeY+bounds.minY, iZ*cellSizeZ+bounds.minZ,
                                        (iX+1)*cellSizeX+bounds.minX, (iY+1)*cellSizeY+bounds.minY, (iZ+1)*cellSizeZ+bounds.minZ };
                grid[iX+iY*cellsX+iZ*cellsX*cellsY] = Cell(cellBounds);
                dimIndex[iX+iY*cellsX+iZ*cellsX*cellsY][0] = iX;
                dimIndex[iX+iY*cellsX+iZ*cellsX*cellsY][1] = iY;
                dimIndex[iX+iY*cellsX+iZ*cellsX*cellsY][2] = iZ;
            }
#endif // 3D
        }
    }
}

void Domain::getNeighborCells(const int &iCell, int *neighborCell){
    // recover iX, iY and iZ
    int iX = dimIndex[iCell][0];
#if DIM >= 2    
    int iY = dimIndex[iCell][1];
#endif // 2D
#if DIM == 3
    int iZ = dimIndex[iCell][2];
#endif // 3D

    int iNeighbor = 0;
    for(int k=iX-1; k<=iX+1; ++k){
#if DIM == 1 
    /// TODO:
#endif // 1D    
#if DIM == 2
    for(int l=iY-1; l<=iY+1; ++l){
            if (k < 0 || k >= cellsX){
                neighborCell[iNeighbor] = -1; // -1: no cell or ghost
            } else if (l < 0 || l >= cellsY) {
                neighborCell[iNeighbor] = -1;
            } else {
                neighborCell[iNeighbor] = k+l*cellsX;
            }
            ++iNeighbor;
#endif // 2D
#if DIM == 3
            for(int m=iZ-1; m<=iZ+1; ++m){
                if (k < 0 || k >= cellsX){
                    neighborCell[iNeighbor] = -1; // -1: no cell or ghost
                } else if (l < 0 || l >= cellsY) {
                    neighborCell[iNeighbor] = -1;
                } else if (m < 0 || m >= cellsZ) {
                    neighborCell[iNeighbor] = -1;
                } else {
                    neighborCell[iNeighbor] = k+l*cellsX+m*cellsX*cellsY;
                }
                ++iNeighbor;
            }
#endif // 3D
        }
    }
}