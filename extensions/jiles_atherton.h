#ifndef _jiles_atherton_h_
#define _jiles_atherton_h_

#include <vector>
#include "cppad/cppad.hpp"

using std::vector;

namespace jiles_atherton {
    template <typename T> vector<T> integrate_RK4_H(vector<double> H, vector<double> t, vector<T> params);
    vector<double> integrate_RK4_H(vector<double> H, vector<double> t, vector<double> params);

    template <typename T> vector<T> integrate_RK4_B(vector<double> B, vector<double> t, vector<T> params);
    vector<double> integrate_RK4_B(vector<double> B, vector<double> t, vector<double> params);

    vector<double> get_gradient(vector<double> H, vector<double> t, vector<double> M, vector<double> params);
    double get_cost(vector<double> H, vector<double> t, vector<double> M, vector<double> params);

    class JilesAthertonModel {
    public:

        JilesAthertonModel();
        JilesAthertonModel(vector<double> H, vector<double> M, vector<double> t);
        ~JilesAthertonModel();

        void calculate(vector<double> params);
        double get_cost(vector<double> params);
        vector<double> get_gradient(vector<double> params);
        //vector< vector<double> > get_hessian();

    private:
        double cost_passive;
        vector<double> H;
        vector<double> M;
        vector<double> t;
        CppAD::ADFun<double> f;
    };
}

#endif