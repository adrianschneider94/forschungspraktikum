import numpy as np


def langevin(x):
    if abs(x) < 1e-3:
        return x / 3.0
    elif abs(x) > 1e8:
        return np.sign(x)
    else:
        return 1 / np.tanh(x) - 1 / x


def grad_langevin(x):
    if abs(x) < 1e-3:
        return 1/3.0 - 1/15.0 * x**2
    elif abs(x) > 1e2:
        return 1/x**2
    elif abs(x) > 1e8:
        return 0.0
    else:
        return 1/x**2 - 1/np.sinh(x)**2
