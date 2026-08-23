
'''Методы для работы с виртуальной температурой'''

from .vec import vec, Number
from .consts import CtoK, M, R_air


@vec()
def e_p_t(e: Number, p: Number, t: Number) -> Number:
    '''
    Расчет виртуальной температуры по парциальному давлению, общему давлению и температуре
    Args:
        e (Number): Парциальное давление водяного пара [Па]
        p (Number): Атмосферное давление [Па]
        t (Number): Температура воздуха [°C]
    Returns:
        Number: Виртуальная температура [K]
    '''
    return (t + CtoK) / (1 - (1 - M) * e / p)


@vec()
def p_rho(p: Number, rho: Number) -> Number:
    '''
    Расчет виртуальной температуры по давлению и плотности
    Args:
        p (Number): Атмосферное давление [Па]
        rho (Number): Плотность влажного воздуха [кг/м³]
    Returns:
        Number: Виртуальная температура [K]
    '''
    return p / (R_air * rho)
