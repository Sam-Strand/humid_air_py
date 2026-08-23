'''Методы для работы с температурой'''

from .vec import vec, Number
from .consts import CtoK, M, k1, k2, C_p_air, C_p_h2o
from math import log


@vec()
def L_sat(L_sat: Number) -> Number:
    '''
    Расчет температуры по удельной теплоте образования пара
    Args:
        L_sat (Number): Удельная теплота образования пара [Дж/кг]
    Returns:
        Number: Температура [°C]
    '''
    return (k1 - L_sat) / k2


@vec()
def e_evap_p_vT(e_evap: Number, p: Number, vT: Number) -> Number:
    '''
    Расчет температуры по парциальному давлению, общему давлению и виртуальной температуре
    Args:
        e_evap (Number): Парциальное давление водяного пара [Па]
        p (Number): Атмосферное давление [Па]
        vT (Number): Виртуальная температура [K]
    Returns:
        Number: Температура воздуха [°C]
    '''
    T = vT * (1 - (1 - M) * e_evap / p)
    return T - CtoK


@vec()
def E_sat(E_sat: Number) -> Number:
    '''
    Расчет температуры по давлению насыщенного пара (обратная формула Магнуса)
    Args:
        E_sat (Number): Давление насыщенного пара [Па]
    Returns:
        Number: Температура воздуха [°C]
    '''
    a = 611.2
    b = 17.62
    c = 243.12
    numerator = log(E_sat / a)
    if numerator < 0:
        b = 22.46
        c = 272.62
        numerator = log(E_sat / a)
    return c * numerator / (b - numerator)


@vec()
def i_d(i: Number, d: Number) -> Number:
    '''
    Расчет температуры по энтальпии и влагосодержанию
    Args:
        i (Number): Энтальпия влажного воздуха [Дж/кг]
        d (Number): Влагосодержание [доля]
    Returns:
        Number: Температура воздуха [°C]
    '''
    return (i - k1 * d) / (C_p_air + d * (C_p_h2o - k2))
