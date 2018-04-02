import numpy as np


def tanh(x):
    return 1.0 - 2.0/(np.exp(2.0*x)+1)

def sinh(x):
    return 0.5*(np.exp(x) - np.exp(-1*x))


def langevin(x):
    if np.abs(x) < 1e-3:
        return x / 3.0
    elif np.abs(x) > 1e8:
        return np.sign(x)
    else:
        return 1 / tanh(x) - 1 / x


def grad_langevin(x):
    if np.abs(x) < 1e-3:
        return 1/3.0 - 1/15.0 * x**2
    elif np.abs(x) > 1e2:
        return 1/x**2
    elif np.abs(x) > 1e8:
        return 0.0
    else:
        return 1/x**2 - 1/sinh(x)**2


def langevin_(x):
    return 1 / tanh(x) - 1 / x


def grad_langevin_(x):
    return 1/x**2 - 1/sinh(x)**2