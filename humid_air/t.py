'''Методы для работы с температурой'''

from .vec import vec, ArrayLike
from .consts import CtoK, M, k1, k2, C_p_air, C_p_h2o
from math import log


@vec(1)
def L(L: ArrayLike):
    '''
    Расчет температуры по удельной теплоте образования пара
    Args:
        L (ArrayLike): Удельная теплота образования пара [Дж/кг]
    Returns:
        ArrayLike: Температура [°C]
    '''
    return (k1 + L) / k2


@vec(3)
def e_evap_p_vT(e_evap: ArrayLike, p: ArrayLike, vT: ArrayLike):
    '''
    Расчет температуры по парциальному давлению, общему давлению и виртуальной температуре
    Args:
        e_evap (ArrayLike): Парциальное давление водяного пара [Па]
        p (ArrayLike): Атмосферное давление [Па]
        vT (ArrayLike): Виртуальная температура [K]
    Returns:
        ArrayLike: Температура воздуха [°C]
    '''
    return (vT / (1 + M * e_evap / p)) - CtoK


@vec(1)
def E_sat(E_sat: ArrayLike):
    '''
    Расчет температуры по давлению насыщенного пара (обратная формула Магнуса)
    Args:
        E_sat (ArrayLike): Давление насыщенного пара [Па]
    Returns:
        ArrayLike: Температура воздуха [°C]
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


@vec(2)
def i_d(i: ArrayLike, d: ArrayLike):
    '''
    Расчет температуры по энтальпии и влагосодержанию
    Args:
        i (ArrayLike): Энтальпия влажного воздуха [Дж/кг]
        d (ArrayLike): Влагосодержание [доля]
    Returns:
        ArrayLike: Температура воздуха [°C]
    '''
    return (i - k1 * d) / (C_p_air + d * (C_p_h2o - k2))
