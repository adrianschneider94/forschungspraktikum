import numpy as np
from forschungspraktikum.functions import langevin, grad_langevin
from scipy.constants import mu_0


def dm_dhe(self, H, M, H_, M_):
    delta = np.sign(H - H_)
    he = H + self.alpha * M
    man = self.m_sat * langevin(he / self.a)
    mirr = (M - self.c * man) / (1.0 - self.c)
    d_mirr_d_he = (man - mirr) / (self.k * delta)
    d_man_d_he = self.m_sat / self.a * grad_langevin(he / self.a)
    return (1 - self.c) * d_mirr_d_he + self.c * d_man_d_he


def dm_dh(self, H, M, H_, M_):
    d_md_he = self.dm_dhe(H, M, H_, M_)
    return d_md_he / (1 - self.alpha * d_md_he)


def dm_db(self, B, M, B_, M_):
    h = B / mu_0 - M
    h_ = B_ / mu_0 - M_
    d_md_he = self.dm_dhe(h, M, h_, M_)
    return 1 / mu_0 * d_md_he / (1 + (1 - self.alpha) * d_md_he)