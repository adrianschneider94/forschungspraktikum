import numpy as np
from forschungspraktikum.jiles_atherton import JilesAthertonModel
from scipy.constants import mu_0
from algopy import CGraph, Function

# Beispielparameter des Jiles-Atherton-Modells
alpha = 0.0021
a = 110.5
k = 30.0
c = 0.4
Msat = 1.35e5

# Eingangsgrößen
r = 2.0e-2 # m, Radius
i_hat = 20.0 # A, Strom
f = 1000.0 # Hz, Frequenz
n = 3 # Anzahl Perioden
n_p = 512# Datenpunkte pro Periode

t = np.arange(n * n_p)/(n_p*f) # Zeitvektor
current = i_hat*(np.sin(2*np.pi*f*t)+0.7*np.sin(6*np.pi*f*t + 1)) # Stromvorgabe
H = current/(2*np.pi*r) # Resultierende Feldvorgabe

graph = CGraph()
graph.trace_on()
x = Function([alpha, a, k, c, Msat])

# Parametervektor
p = {'alpha': x[0],
     'a': x[1],
     'k': x[2],
     'c': x[3],
     'm_sat': x[4]}

model = JilesAthertonModel.from_dict(p)
M = model.integrate_rk4(t, H)

H = H[::2]
t = t[::2]
B = mu_0*(H + M)
dB_dt = np.zeros(np.size(B))
new = np.append([0.0], (B[1:] - B[0:-1]) / (t[1:] - t[0:-1]))

P = np.sum(0.5*H*new)

graph.trace_off()
graph.independentFunctionList = [x]
graph.dependentFunctionList = [P]

a = graph.gradient([alpha, a, k, c, Msat])
print(a)