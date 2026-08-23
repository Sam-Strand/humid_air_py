'''Методы для работы с парциальным давлением водяного пара'''

from .vec import vec, Number
from .consts import M, CtoK


@vec()
def d_p(d: Number, p: Number) -> Number:
    '''
    Расчет парциального давления по влагосодержанию и общему давлению
    Args:
        d (Number): Влагосодержание [доля]
        p (Number): Атмосферное давление [Па]
    Returns:
        Number: Парциальное давление водяного пара [Па]
    '''
    return (d * p) / (M + d)


@vec()
def E_sat_h(E_sat: Number, h: Number) -> Number:
    '''
    Расчет парциального давления по давлению насыщения и влажности
    Args:
        E_sat (Number): Давление насыщенного пара [Па]
        h (Number): Относительная влажность [доля]
    Returns:
        Number: Парциальное давление водяного пара [Па]
    '''
    return h * E_sat


@vec()
def p_t_vT(p: Number, t: Number, vT: Number) -> Number:
    '''
    Расчет парциального давления по общему давлению, температуре и виртуальной температуре
    Args:
        p (Number): Атмосферное давление [Па]
        t (Number): Температура воздуха [°C]
        vT (Number): Виртуальная температура [K]
    Returns:
        Number: Парциальное давление водяного пара [Па]
    '''
    T = t + CtoK
    return (p / (1 - M)) * (1 - T / vT)
