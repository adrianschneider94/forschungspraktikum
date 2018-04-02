#ifndef _functions_h_
#define _functions_h_

#include "cppad/cppad.hpp"

using CppAD::AD;

namespace functions{
    template<typename T> T langevin(T x);
    template<typename T> T grad_langevin(T x);
    template<typename T> T sign(T x);

    double grad_langevin(double x);
    double langevin(double x);
    double sign(double x);

    AD<double> grad_langevin(AD<double> x);
    AD<double> langevin(AD<double> x);
    AD<double> sign(AD<double> x);
}

#endif