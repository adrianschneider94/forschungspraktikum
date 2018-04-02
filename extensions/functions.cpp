//
// Created by Adrian Schneider on 10.03.18.
//

#include "cppad/cppad.hpp"

using CppAD::AD;

namespace functions {

    template<typename T>
    T sign(T x) {
        if (x >= 0) {
            return 1.0;
        } else {
            return -1.0;
        }
    }

    template<typename T>
    T langevin(T x) {
        if (CppAD::abs(x) < 1.0e-8) {
            return x/3.0;
        } else if (CppAD::abs(x) > 1.0e8){
            return CppAD::sign(x);
        }
        else {
            return 1 / CppAD::tanh(x) - 1 / x;
        }
    }

    template<typename T>
    T grad_langevin(T x) {
        if (CppAD::abs(x) < 1.0e-8){
            return 1.0/3.0 - 1.0/15.0 * CppAD::pow(x, 2);
        } else if (CppAD::abs(x) > 1.0e8){
            return CppAD::pow(x, -2.0);
        } else if (CppAD::abs(x) > 1.0e8) {
            return 0.0;
        } else {
            return CppAD::pow(x, -2.0) - CppAD::pow(CppAD::sinh(x), -2.0);
        }
    }

    double sign(double x) {
        return  sign<double>(x);
    }

    double langevin(double x) {
        return langevin<double>(x);
    }

    double grad_langevin(double x) {
        return grad_langevin<double>(x);
    }

    AD<double> sign(AD<double> x) {
        return  sign< AD<double> >(x);
    }

    AD<double> langevin(AD<double> x) {
        return langevin< AD<double> >(x);
    }

    AD<double> grad_langevin(AD<double> x) {
        return grad_langevin< AD<double> >(x);
    }
}

