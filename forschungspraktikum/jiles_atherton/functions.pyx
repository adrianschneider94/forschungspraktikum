# distutils: language = c++
# distutils: include_dirs = ./extensions
# distutils: sources = ./extensions/functions.cpp ./extensions/jiles_atherton.cpp

import numpy as np
cimport numpy as np

cimport functions_c

def langevin(x):
    return functions_c.langevin(x)

def grad_langevin(x):
    return functions_c.grad_langevin(x)

def jiles_atherton_H(x, y, z):
    return functions_c.integrate_RK4_H(x, y, z)

def jiles_atherton_B(x, y, z):
    return functions_c.integrate_RK4_B(x, y, z)

def get_gradient(H, t, M, params):
    return functions_c.get_gradient(H, t, M, params)

def get_cost(H, t, M, params):
    return functions_c.get_cost(H, t, M, params)

cdef class JilesAthertonModel:
    cdef functions_c.JilesAthertonModel c_JilesAthertonModel

    def __cinit__(self, H, M, t):
        self.c_JilesAthertonModel = functions_c.JilesAthertonModel(H, M, t)

    def calculate(self, params):
        self.c_JilesAthertonModel.calculate(params)

    def get_cost(self, params):
        return self.c_JilesAthertonModel.get_cost(params)

    def get_gradient(self, params):
        return self.c_JilesAthertonModel.get_gradient(params)

