
'''Методы для работы с виртуальной температурой'''

from .vec import vec, ArrayLike
from .consts import CtoK, M, R_air


@vec(3)
def e_p_t(e: ArrayLike, p: ArrayLike, t: ArrayLike):
    '''
    Расчет виртуальной температуры по парциальному давлению, общему давлению и температуре
    Args:
        e (ArrayLike): Парциальное давление водяного пара [Па]
        p (ArrayLike): Атмосферное давление [Па]
        t (ArrayLike): Температура воздуха [°C]
    Returns:
        ArrayLike: Виртуальная температура [K]
    '''
    return (t + CtoK) / (1 - (1 - M) * e / p)


@vec(2)
def p_rho(p: ArrayLike, rho: ArrayLike):
    '''
    Расчет виртуальной температуры по давлению и плотности
    Args:
        p (ArrayLike): Атмосферное давление [Па]
        rho (ArrayLike): Плотность влажного воздуха [кг/м³]
    Returns:
        ArrayLike: Виртуальная температура [K]
    '''
    return p / (R_air * rho)
