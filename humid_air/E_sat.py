'''Методы для работы с давлением насыщенного пара'''

from .vec import vec, ArrayLike
from math import exp


@vec(1)
def t(t: ArrayLike):
    '''
    Расчет давления насыщенного пара по температуре (формула Магнуса)
    Args:
        t (ArrayLike): Температура воздуха [°C]
    Returns:
        ArrayLike: Давление насыщенного пара [Па]
    '''
    a = 611.2
    b = 17.62 if t > 0 else 22.46
    c = 243.12 if t > 0 else 272.62
    return a * exp(b * t / (c + t))


@vec(2)
def e_evap_h(e_evap: ArrayLike, h: ArrayLike):
    '''
    Расчет давления насыщенного пара по парциальному давлению и влажности
    Args:
        e_evape (ArrayLike): Парциальное давление водяного пара [Па]
        h (ArrayLike): Относительная влажность [доля]
    Returns:
        ArrayLike: Давление насыщенного пара [Па]
    '''
    return e_evap / h
