# distutils: language = c++
# distutils: include_dirs = ./extensions
# distutils: sources = ./extensions/jiles_atherton.cpp

from libcpp.vector cimport vector

cdef extern from "functions.h" namespace "functions":
    cdef double langevin(double)
    cdef double grad_langevin(double)

cdef extern from "jiles_atherton.h" namespace "jiles_atherton":
    cdef vector[double] integrate_RK4_H(vector[double], vector[double], vector[double])
    cdef vector[double] integrate_RK4_B(vector[double], vector[double], vector[double])

    cdef vector[double] get_gradient(vector[double], vector[double], vector[double], vector[double])
    cdef double get_cost(vector[double], vector[double], vector[double], vector[double])

    cdef cppclass JilesAthertonModel:
        JilesAthertonModel() except +
        JilesAthertonModel(vector[double], vector[double], vector[double]) except +
        void calculate(vector[double])
        double get_cost(vector[double])
        vector[double] get_gradient(vector[double])

