'''Методы для работы с парциальным давлением водяного пара'''

from .vec import vec, ArrayLike
from .consts import M, CtoK


@vec(2)
def d_p(d: ArrayLike, p: ArrayLike):
    '''
    Расчет парциального давления по влагосодержанию и общему давлению
    Args:
        d (ArrayLike): Влагосодержание [доля]
        p (ArrayLike): Атмосферное давление [Па]
    Returns:
        ArrayLike: Парциальное давление водяного пара [Па]
    '''
    return (d * p) / (M + d)


@vec(2)
def E_sat_h(E_sat: ArrayLike, h: ArrayLike):
    '''
    Расчет парциального давления по давлению насыщения и влажности
    Args:
        E_sat (ArrayLike): Давление насыщенного пара [Па]
        h (ArrayLike): Относительная влажность [доля]
    Returns:
        ArrayLike: Парциальное давление водяного пара [Па]
    '''
    return h * E_sat


@vec(3)
def p_t_vT(p: ArrayLike, t: ArrayLike, vT: ArrayLike):
    '''
    Расчет парциального давления по общему давлению, температуре и виртуальной температуре
    Args:
        p (ArrayLike): Атмосферное давление [Па]
        t (ArrayLike): Температура воздуха [°C]
        vT (ArrayLike): Виртуальная температура [K]
    Returns:
        ArrayLike: Парциальное давление водяного пара [Па]
    '''
    T = t + CtoK
    return (p / (1 - M)) * (1 - T / vT)
