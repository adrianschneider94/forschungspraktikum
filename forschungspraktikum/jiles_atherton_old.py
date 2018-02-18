import numpy as np
from scipy.constants import mu_0




class JilesAthertonModel(object):
    def __init__(self, alpha=0.0, m_sat=0.0, a=0.0, c=0.0, k=0.0):
        super(JilesAthertonModel, self).__init__()
        self.alpha = alpha
        self.a = a
        self.c = c
        self.k = k
        self.m_sat = m_sat

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