'''Методы для работы с влагосодержанием'''

from .vec import vec, ArrayLike
from .consts import C_p_air, M, C_p_h2o, k1, k2


@vec(2)
def e_evap_p(e_evap: ArrayLike, p: ArrayLike):
    '''
    Расчет влагосодержания по парциальному давлению и общему давлению
    Args:
        e_evap (ArrayLike): Парциальное давление водяного пара [Па]
        p (ArrayLike): Атмосферное давление [Па]
    Returns:
        ArrayLike: Влагосодержание [доля]
    '''
    return M * e_evap / (p - e_evap)


@vec(2)
def L_sat_l_evap(L_sat: ArrayLike, l_evap: ArrayLike):
    '''
    Расчет влагосодержания по удельной и общей теплоте парообразования
    Args:
        L_sat (ArrayLike): Удельная теплота образования пара [Дж/кг]
        l_evap (ArrayLike): Удельная теплота парообразования [Дж/кг]
    Returns:
        ArrayLike: Влагосодержание [доля]
    '''
    return l_evap / L_sat


@vec(1)
def c(c: ArrayLike):
    '''
    Расчет влагосодержания по удельной теплоемкости
    Args:
        c (ArrayLike): Удельная теплоемкость влажного воздуха [Дж/(кг·K)]
    Returns:
        ArrayLike: Влагосодержание [доля]
    '''
    return (c - C_p_air) / (C_p_h2o - C_p_air)


@vec(2)
def i_t(i: ArrayLike, t: ArrayLike):
    '''
    Расчет влагосодержания по энтальпии и температуре
    Args:
        i (ArrayLike): Энтальпия влажного воздуха [Дж/кг]
        t (ArrayLike): Температура воздуха [°C]
    Returns:
        ArrayLike: Влагосодержание [доля]
    '''
    return (i - C_p_air * t) / (k1 + (C_p_h2o - k2) * t)
