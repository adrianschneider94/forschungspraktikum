//
// Created by Adrian Schneider on 10.03.18.
//
#include <iostream>
#include <vector>
#include <cmath>
#include "functions.h"
#include "jiles_atherton.h"

double mu_0_ = 1.25663706e-6;

int main() {
    vector<double> params(5);
    params[0] = 0.000;
    params[1] = 29.5;
    params[2] = 387000;
    params[3] = 30.0;
    params[4] = 0.3;

    size_t n = 8192 * 4;
    size_t i;
    vector<double> B(n);
    vector<double> time(n);
    vector<double> H(n);
    for (i = 0; i < n; i++){
        time[i] = (double) i/n;
        B[i] = 0.4*std::cos(8 * M_PI * time[i])/mu_0_;
    }
    H = jiles_atherton::integrate_RK4_B(B, time, params);
    size_t j;
    for (j=0; j < H.size(); j++) {
        std::cout << H[i] << std::endl;
    }
    return 0;
}