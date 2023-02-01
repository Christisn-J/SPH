//
// Created by Christian Jetter on 23.01.23
//

#include "../include/kernel.h"

double Kernel::cubicSpline(const double &r, const double &h) {

#if DIM == 1 //TODO 
    const double sigma = nan;
#endif // 1D

#if DIM == 2
    const double sigma = 40./(7.*M_PI*pow(h,DIM));
#endif // 2D

#if DIM == 3
    const double sigma = 8./(M_PI*pow(h,DIM));
#endif // 3D

    const double q = r/h;
    if (0. <= q && q < 0.5){
        return sigma*(6.*pow(q,3)-6.*pow(q,2)+1);
    } else if (0.5 <= q && q <= 1.){
        return sigma*2.*pow(1.-q, 3.);
    } else {
        return 0.;
    }
}