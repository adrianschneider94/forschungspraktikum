import numpy as np
from scipy.misc import factorial
from algopy import CGraph, Function


def cos_series(x):
    res = 0
    for n in range(10):
        res += (-1)**n*x**(2*n)/factorial(2*n)
    return res


graph = CGraph()
graph.trace_on()
x = Function(2.0)

res = cos_series(x)

graph.trace_off()
graph.independentFunctionList = [x]
graph.dependentFunctionList = [res]

a = graph.gradient([3/4*np.pi])
print(a)