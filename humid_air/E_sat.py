'''Методы для работы с давлением насыщенного пара'''

from .vec import vec, Number
from numpy import exp


@vec()
def t(t: Number) -> Number:
    '''
    Расчет давления насыщенного пара по температуре (формула Магнуса)
    Args:
        t (Number): Температура воздуха [°C]
    Returns:
        Number: Давление насыщенного пара [Па]
    '''
    a = 611.2
    b = 17.62 if t > 0 else 22.46
    c = 243.12 if t > 0 else 272.62
    return a * exp(b * t / (c + t))


@vec()
def e_evap_h(e_evap: Number, h: Number) -> Number:
    '''
    Расчет давления насыщенного пара по парциальному давлению и влажности
    Args:
        e_evape (Number): Парциальное давление водяного пара [Па]
        h (Number): Относительная влажность [доля]
    Returns:
        Number: Давление насыщенного пара [Па]
    '''
    return e_evap / h
