import numpy as np
from scipy.misc import factorial
from algopy import CGraph, Function


class AlgopyTest(object):
    def __init__(self, a):
        super(AlgopyTest, self).__init__()
        self.a = a

    def cos_series(self, x):
        res = 0
        for n in range(10):
            res += (-1) ** n * x ** (2 * n) / factorial(2 * n)
        return self.a*res


graph = CGraph()
graph.trace_on()
x = Function([2.0, 1.0])
test = AlgopyTest(x[1])
res = test.cos_series(x[0])

graph.trace_off()
graph.independentFunctionList = [x]
graph.dependentFunctionList = [res]

a = graph.gradient([3/4*np.pi, 4.0])
print(a)