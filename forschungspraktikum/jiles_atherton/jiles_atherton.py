import numpy as np
from forschungspraktikum.jiles_atherton.functions_ import langevin_ as langevin
from forschungspraktikum.jiles_atherton.functions_ import grad_langevin_ as grad_langevin


class JilesAthertonModel(object):

    def __init__(self, a, alpha, k, m_sat, c):
        super(JilesAthertonModel, self).__init__()
        self.a = a
        self.alpha = alpha
        self.k = k
        self.m_sat = m_sat
        self.c = c
        self.params = a, alpha, k, m_sat, c

    @classmethod
    def from_dict(cls, params):
        if type(params) is not dict:
            raise TypeError("Das Argument p muss ein Dictionary sein.")

        keys = ['a', 'alpha', 'm_sat', 'k', 'c']
        if False in map(lambda key: key in params, keys):
            raise ValueError("Der Parameter p enthält nicht alle notwendigen Schlüssel.")
        return cls(params['a'], params['alpha'], params['k'], params['m_sat'], params['c'])

    def anhysteric_magnetization(self, he):
        return self.m_sat * langevin(he / self.a)

    def d_anhysteric_magnetization_wrt_effective_magnetic_field(self, he):
        return self.m_sat / self.a * grad_langevin(he / self.a)

    def d_irreversible_magnetization_wrt_effective_magnetic_field(self, man, mirr, dh_dt):
        return (man - mirr) / (self.k * np.sign(dh_dt))

    def d_magnetization_wrt_effective_magnetic_field(self, m, h, dh_dt):
        he = h + self.alpha * m
        m_an = self.anhysteric_magnetization(he)
        m_irr = (m - self.c * m_an) / (1.0 - self.c)
        dmirr_dhe = self.d_irreversible_magnetization_wrt_effective_magnetic_field(m_an, m_irr, dh_dt)
        dman_dhe = self.d_anhysteric_magnetization_wrt_effective_magnetic_field(he)

        return (1.0 - self.c) * dmirr_dhe + self.c * dman_dhe

    def d_magnetization_wrt_magnetic_field(self, m, h, dh_dt):
        dm_dhe = self.d_magnetization_wrt_effective_magnetic_field(m, h, dh_dt)
        return dm_dhe / (1.0 - self.alpha * dm_dhe)

    def integrate_explicit_euler(self, t, H):
        M = np.zeros(np.size(t))

        # Zur einfacheren Notation: Funktion umbenennen
        JA = self.d_magnetization_wrt_magnetic_field

        for i in range(1, np.size(t)):
            dH_dt = (H[i - 1] - H[i - 2]) / (t[i - 1] - t[i - 2])
            M[i] = M[i - 1] + JA(M[i - 1], H[i - 1], dH_dt) * dH_dt * (t[i] - t[i - 1])

        return M

    def integrate_rk4(self, t, H):
        # Vektor für Magnetisierung initialisieren, doppelte Schrittweite wie H-Feld!
        M = [0.0, 0.0]

        # Zur Vereinfachung: dH/dt als Vektor aufstellen
        dH_dt = np.zeros(np.size(H))
        dH_dt[1:] = (H[1:] - H[0:-1]) / (t[1:] - t[0:-1])

        # Zur einfacheren Notation: Funktion umbenennen
        JA = self.d_magnetization_wrt_magnetic_field

        # Runge-Kutta-Iteration
        for n in range(2, int(np.size(H)/2)):
            h = t[2 * n] - t[2 * n - 2]

            k1 = JA(M[n - 1], H[2 * n - 2], dH_dt[2 * n - 2]) * dH_dt[2 * n - 2]
            k2 = JA(M[n - 1] + h * k1 / 2.0, H[2 * n - 1], dH_dt[2 * n - 1]) * dH_dt[2 * n - 1]
            k3 = JA(M[n - 1] + h * k2 / 2.0, H[2 * n - 1], dH_dt[2 * n - 1]) * dH_dt[2 * n - 1]
            k4 = JA(M[n - 1] + h * k3, H[2 * n], dH_dt[2 * n]) * dH_dt[2 * n]

            x = M[n - 1] + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
            M.append(x)

        return M
