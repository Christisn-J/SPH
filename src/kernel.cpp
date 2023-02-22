//
// Created by Christian Jetter on 23.01.23
//

#include "../include/kernel.h"

double KernelR::cubicSpline(const double &r, const double &h) {
#if DIM == 1 //TODO 
    const double sigma = nan;
#endif // 1D

#if DIM == 2
    const double sigma = 40./(7.*M_PI);
#endif // 2D

#if DIM == 3
    const double sigma = 8./(M_PI);
#endif // 3D
    
    double o = sigma/pow(h,DIM);

    const double q = r/h;
    if (0. <= q && q < 0.5){
        return o*(6.*pow(q,3)-6.*pow(q,2)+1);
    } else if (0.5 <= q && q <= 1.){
        return o*2.*pow(1.-q, 3.);
    } else {
#if DEBUG_LVL == NOT_IN_USE
        Logger(DEBUG) << "Kernal set to 0.0";
#endif
        return 0.;
    }
}

double NablaKernelR::cubicSpline(const double &r, const double &h) {
#if DIM == 1 //TODO 
    const double sigma = nan;
#endif // 1D

#if DIM == 2
    const double sigma = 40./(7.*M_PI);
#endif // 2D

#if DIM == 3
    const double sigma = 8./(M_PI);
#endif // 3D

    double o = 6*sigma/pow(h,DIM+1);
    const double q = r/h;
    if (0. <= q && q < 0.5){
        return o*(3.*pow(q,2.)-2.*q);
    } else if (0.5 <= q && q <= 1.){
        return o*2.*-pow(1.-q, 2.);
    } else {
#if DEBUG_LVL == NOT_IN_UESED
        Logger(DEBUG) << "dKernal set to 0.0";
#endif
        return 0.;
    }
}