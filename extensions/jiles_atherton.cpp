//
// Created by Adrian Schneider on 10.03.18.
//
#include <vector>
#include <cmath>
#include <iostream>
#include "functions.h"
#include "jiles_atherton.h"
#include "cppad/cppad.hpp"

using std::vector;
double mu_0 = 1.25663706e-6;

namespace jiles_atherton {
    template<typename T>
    T anhysteric_magnetization(T H_e, T a, T M_sat) {
        T M_an;
        M_an = M_sat * functions::langevin(H_e / a);
        return M_an;
    }

    template<typename T>
    T d_anhysteric_magnetization_wrt_effective_magnetic_field(T H_e, T a, T M_sat) {
        T d_Man_d_He;
        d_Man_d_He = M_sat / a * functions::grad_langevin(H_e / a);
        return d_Man_d_He;
    }

    template<typename T>
    T d_irreversible_magnetization_wrt_effective_magnetic_field(T M_an, T M_irr, T dH_dt, T k) {
        T d_Mirr_d_He;
        d_Mirr_d_He = (M_an - M_irr) / (k * functions::sign(dH_dt));
        return d_Mirr_d_He;
    }

    template<typename T>
    T d_magnetization_wrt_effective_magnetic_field(T M, T H, T dH_dt, T alpha, T a, T M_sat, T k, T c) {
        T H_e, M_an, M_irr, dMirr_dHe, dMan_dHe, dM_dHe;
        H_e = H + alpha * M;
        M_an = anhysteric_magnetization<T>(H_e, a, M_sat);
        M_irr = (M - c * M_an) / (1.0 - c);
        dMirr_dHe = d_irreversible_magnetization_wrt_effective_magnetic_field<T>(M_an, M_irr, dH_dt, k);
        dMan_dHe = d_anhysteric_magnetization_wrt_effective_magnetic_field<T>(H_e, a, M_sat);

        dM_dHe = (1.0 - c) * dMirr_dHe + c * dMan_dHe;
        return dM_dHe;
    }

    template<typename T>
    T d_magnetization_wrt_magnetic_field(T M, T H, T dH, vector<T> params) {
        T dM_dHe, alpha, a, M_sat, k, c;

        alpha = params[0];
        a = params[1];
        M_sat = params[2];
        k = params[3];
        c = params[4];

        dM_dHe = d_magnetization_wrt_effective_magnetic_field<T>(M, H, dH, alpha, a, M_sat, k, c);
        return 1.0 / (1.0/dM_dHe - alpha);
    }

    template<typename T>
    T d_magnetic_field_wrt_magnetic_flux(T H, T B, T dB, vector<T> params) {
        T dM_dHe, alpha, a, M_sat, k, c, M, dM_dB;

        alpha = params[0];
        a = params[1];
        M_sat = params[2];
        k = params[3];
        c = params[4];

        M = B/mu_0 - H;

        dM_dHe = d_magnetization_wrt_effective_magnetic_field(M, H, dB, alpha, a, M_sat, k, c);
        dM_dB = 1/mu_0 * dM_dHe/(1.0 + (1.0-alpha)*dM_dHe);

        return 1/mu_0*((1-alpha*dM_dHe)/(1 + (1-alpha)*dM_dHe));
    }

    template<typename T>
    vector<T> integrate_RK4_H(vector<double> H, vector<double> t, vector<T> params) {
        // Construct the derivative dH
        vector<double> dH(H.size());
        size_t i;
        for (i = 1; i < H.size(); i++) {
            dH[i] = H[i] - H[i - 1];
        }

        // Initialize the vector for the magnetization
        size_t M_size = H.size() / 2;
        vector<T> M(M_size);

        // Runge-Kutta iteration
        size_t n;
        T h, k1, k2, k3, k4;
        for (n = 2; n < M_size; n++) {
            h = H[2*n] - H[2*n - 2];

            k1 = d_magnetization_wrt_magnetic_field<T>(M[n - 1], H[2 * n - 2], dH[2 * n - 2], params) * h;
            k2 = d_magnetization_wrt_magnetic_field<T>(M[n - 1] + k1 / 2.0, H[2 * n - 1], dH[2 * n - 1], params) * h;
            k3 = d_magnetization_wrt_magnetic_field<T>(M[n - 1] + k2 / 2.0, H[2 * n - 1], dH[2 * n - 1], params) * h;
            k4 = d_magnetization_wrt_magnetic_field<T>(M[n - 1] + k3, H[2 * n], dH[2 * n], params) * h;

            M[n] = M[n - 1] + 1.0 / 6.0 * (k1 + 2.0 * k2 + 2.0 * k3 + k4);

        }

        return M;
    }

    template<typename T>
    vector<T> integrate_RK4_B(vector<double> B, vector<double> t, vector<T> params) {
        // Construct the derivative dH
        vector<double> dB(B.size());
        size_t i;
        for (i = 1; i < B.size(); i++) {
            dB[i] = B[i] - B[i - 1];
        }

        vector<T> H(B.size()/2);

        T alpha, a, M_sat, k, c;

        alpha = params[0];
        a = params[1];
        M_sat = params[2];
        k = params[3];
        c = params[4];

        size_t n;
        T h, k1, k2, k3, k4;
        for (n = 1; n < H.size(); n++) {
            h = B[2*n] - B[2*n - 2];

            k1 = d_magnetic_field_wrt_magnetic_flux<T>(H[n - 1], B[2 * n - 2], dB[2 * n - 2], params) * h;
            k2 = d_magnetic_field_wrt_magnetic_flux<T>(H[n - 1] + k1 / 2.0, B[2 * n - 1], dB[2 * n - 1], params) * h;
            k3 = d_magnetic_field_wrt_magnetic_flux<T>(H[n - 1] + k2 / 2.0, B[2 * n - 1], dB[2 * n - 1], params) * h;
            k4 = d_magnetic_field_wrt_magnetic_flux<T>(H[n - 1] + k3, B[2 * n], dB[2 * n], params) * h;
            std::cout << k1 << " " << k2 << " " << k3 << " " << k4 <<std::endl;

            H[n] = H[n - 1] + 1.0 / 6.0 * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
        }

        return H;
    }


    vector<double> integrate_RK4_H(vector<double> H, vector<double> t, vector<double> params) {
        return integrate_RK4_H<double>(H, t, params);
    }


    vector<double> integrate_RK4_B(vector<double> B, vector<double> t, vector<double> params) {
        return integrate_RK4_B<double>(B, t, params);
    }


    vector<double> get_gradient(vector<double> H, vector<double> t, vector<double> M, vector<double> params) {
        using CppAD::AD;
        vector<AD<double> > a_params(5);
        vector<AD<double> > cost(1);
        size_t i;
        for (i = 0; i < params.size(); i++) {
            a_params[i] = params[i];
        }

        CppAD::Independent(a_params);

        vector<AD<double> > result(H.size() / 2);
        result = integrate_RK4_H<AD<double> >(H, t, a_params);
        cost[0] = 0;
        for (i = 0; i < result.size(); i++) {
            cost[0] += CppAD::pow(result[i] - M[i], 2) / (2 * H.size());
        }

        CppAD::ADFun<double> f;
        f.Dependent(a_params, cost);
        f.optimize();

        vector<double> jac(params.size());
        jac = f.Jacobian(params);

        return jac;
    }

    double get_cost(vector<double> H, vector<double> t, vector<double> M, vector<double> params) {
        double cost;
        vector<double> result(H.size() / 2);
        result = integrate_RK4_H<double>(H, t, params);
        cost = 0;
        size_t i;
        for (i = 0; i < result.size(); i++) {
            cost += CppAD::pow(result[i] - M[i], 2) / (2 * H.size());
        }
        return cost;
    }

    JilesAthertonModel::JilesAthertonModel() {}


    JilesAthertonModel::JilesAthertonModel(vector<double> H, vector<double> M, vector<double> t)
    : H(H), M(M), t(t) {}

    JilesAthertonModel::~JilesAthertonModel() {}

    void JilesAthertonModel::calculate(vector<double> params_) {
        using CppAD::AD;
        vector< AD<double> > a_params(params_.size());
        vector< AD<double> > cost(1);

        size_t i;
        for (i = 0; i < params_.size(); i++) {
            a_params[i] = params_[i];
        }

        CppAD::Independent(a_params);

        vector< AD<double> > result(H.size() / 2);
        result = integrate_RK4_H<AD<double> >(H, t, a_params);

        cost[0] = 0;
        for (i = 0; i < result.size(); i++) {
            cost[0] += CppAD::pow(result[i] - M[i], 2) / (2 * H.size());
        }

        f.Dependent(a_params, cost);
        f.optimize();
    }

    double JilesAthertonModel::get_cost(vector<double> params) {
        return f.Forward(0, params)[0];
    }

    vector<double> JilesAthertonModel::get_gradient(vector<double> params) {
        vector<double> jac(params.size());
        jac = f.Jacobian(params);

        return jac;
    }
}