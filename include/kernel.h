//
// Created by Christian Jetter on 23.01.23
//

#ifndef KERNAL_H
#define KERNAL_H

#include <cmath>
#include "global.h"
#include "Logger.h"

namespace KernelR {
    double cubicSpline(const double &r, const double &h);
}

namespace NablaKernelR {
    double cubicSpline(const double &r, const double &h);
}


#endif // KERNAL_H