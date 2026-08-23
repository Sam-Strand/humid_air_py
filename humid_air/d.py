'''Методы для работы с влагосодержанием'''

from .vec import vec, Number
from .consts import C_p_air, M, C_p_h2o, k1, k2


@vec()
def e_evap_p(e_evap: Number, p: Number) -> Number:
    '''
    Расчет влагосодержания по парциальному давлению и общему давлению
    Args:
        e_evap (Number): Парциальное давление водяного пара [Па]
        p (Number): Атмосферное давление [Па]
    Returns:
        Number: Влагосодержание [доля]
    '''
    return M * e_evap / (p - e_evap)


@vec()
def L_sat_l_evap(L_sat: Number, l_evap: Number) -> Number:
    '''
    Расчет влагосодержания по удельной и общей теплоте парообразования
    Args:
        L_sat (Number): Удельная теплота образования пара [Дж/кг]
        l_evap (Number): Удельная теплота парообразования [Дж/кг]
    Returns:
        Number: Влагосодержание [доля]
    '''
    return l_evap / L_sat


@vec()
def c(c: Number) -> Number:
    '''
    Расчет влагосодержания по удельной теплоемкости
    Args:
        c (Number): Удельная теплоемкость влажного воздуха [Дж/(кг·K)]
    Returns:
        Number: Влагосодержание [доля]
    '''
    return (c - C_p_air) / (C_p_h2o - C_p_air)


@vec()
def i_t(i: Number, t: Number) -> Number:
    '''
    Расчет влагосодержания по энтальпии и температуре
    Args:
        i (Number): Энтальпия влажного воздуха [Дж/кг]
        t (Number): Температура воздуха [°C]
    Returns:
        Number: Влагосодержание [доля]
    '''
    return (i - C_p_air * t) / (k1 + (C_p_h2o - k2) * t)
