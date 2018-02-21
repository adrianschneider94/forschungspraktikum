from numpy import sign
from forschungspraktikum.functions import langevin, grad_langevin

""" Jiles-Atherton-Modell

In diesem Modul werden die Gleichungen des Jiles-Atherton-Modells implementiert. Das Jiles-Atherton-Modell hat fünf 
Parameter, die in diesem Modul jeweils in einen Parametervektor zusammengefasst werden. Es gilt stets:


Notes
_____
Für den Parametervektor gilt:
    p['alpha']  :math:`\alpha`          Interdomänenkopplung
    p['a']      :math:`a`               Domänenwanddichte
    p['m_sat']  :math:`M_\text{sat}`    Sättigungsmagnetisierung
    p['k']      :math:`k`               Pinning-Energie
    p['c']      :math:`c`               Magnetisierungsreversibilität

"""


def anhysteric_magnetization(he, p):
    """

    Parameters
    ----------
    he : float
        Die effektive magnetische Feldstärke H.
    p : dict
        Der Parametervektor p.

    Returns
    -------
    float
        Die anhysterischen Magnetisierung.

    """
    a = p['a']
    m_sat = p['m_sat']

    return m_sat * langevin(he/a)


def d_anhysteric_magnetization_wrt_effective_magnetic_field(he, p):
    """

    Parameters
    ----------
    he : float
        Die effektive magnetische Feldstärke H.
    p : dict
        Der Parametervektor p.

    Returns
    -------
    float
        Das Differential der anhysterischen Magnetisierung nach der effektiven magnetischen Feldstärke H_e, dM_an/dH_e.

    """
    a = p['a']
    m_sat = p['m_sat']

    return m_sat/a * grad_langevin(he/a)


def d_irreversible_magnetization_wrt_effective_magnetic_field(man, mirr, dh_dt, p):
    """

    Parameters
    ----------
    man : float
        Die anhysterische Magnetisierung.
    mirr : float
        Die irreversible Magnetisierung.
    dh_dt : float
        Die zeitliche Ableitung der magnetischen Feldstärke.
    p : dict
        Der Parametervektor p.

    Returns
    -------
    float
        Das Differential der irreversiblen Magnetisierung nach der effektiven magnetischen Feldstärke H_e, dM_irr/dH_e.

    """
    k = p['k']
    return (man - mirr)/(k * sign(dh_dt))


def d_magnetization_wrt_effective_magnetic_field(m, h, dh_dt, p):
    """

    Parameters
    ----------
    h : float
        Die magnetische Feldstärke H.
    m : float
        Die Magnetisierung M.
    dh_dt : float
        Die zeitliche Ableitung der magnetischen Feldstärke.
    p : dict
        Der Parametervektor p.

    Returns
    -------
    float
        Das Differential der Magnetisierung nach der effektiven magnetischen Feldstärke H_e, dM/dH_e.

    """
    alpha = p['alpha']
    c = p['c']

    he = h + alpha * m
    m_an = anhysteric_magnetization(he, p)
    m_irr = (m - c * m_an)/(1.0 - c)
    dmirr_dhe = d_irreversible_magnetization_wrt_effective_magnetic_field(m_an, m_irr, dh_dt, p)
    dman_dhe = d_anhysteric_magnetization_wrt_effective_magnetic_field(he, p)

    return (1.0 - c) * dmirr_dhe + c * dman_dhe


def d_magnetization_wrt_magnetic_field(m, h, dh_dt, p):
    """

    Parameters
    ----------
    h : float
        Die magnetische Feldstärke H.
    m : float
        Die Magnetisierung M.
    dh_dt : float
        Die zeitliche Ableitung der magnetischen Feldstärke.
    p : dict
        Der Parametervektor p.

    Returns
    -------
    float
        Das Differential der Magnetisierung nach der magnetischen Feldstärke H, dM/dH.

    """
    alpha = p['alpha']

    dm_dhe = d_magnetization_wrt_effective_magnetic_field(m, h, dh_dt, p)

    return dm_dhe/(1.0 - alpha * dm_dhe)


JilesAthertonH = d_magnetization_wrt_magnetic_field
