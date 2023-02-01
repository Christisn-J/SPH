//
// Created by Christian Jetter on 23.01.23
//

#include "../include/lib.h"

double distance2D(const double xi,const double yi,const double xj, const double yj){
    return sqrt(pow(xi-xj,2)+pow(yi-yj,2));
}