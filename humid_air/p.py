'''Методы для работы с атмосферным давлением'''

from .vec import vec, Number
from .consts import CtoK, M, R_air


@vec()
def e_t_vT(e: Number, t: Number, vT: Number) -> Number:
    '''
    Расчет общего давления по парциальному давлению, температуре и виртуальной температуре
    Args:
        e (Number): Парциальное давление водяного пара [Па]
        t (Number): Температура воздуха [°C]
        vT (Number): Виртуальная температура [K]
    Returns:
        Number: Атмосферное давление [Па]
    '''
    T = t + CtoK
    return (1 - M) * e * vT / (vT - T)


@vec()
def d_e(d: Number, e: Number) -> Number:
    '''
    Расчет общего давления по влагосодержанию и парциальному давлению
    Args:
        d (Number): Влагосодержание [доля]
        e (Number): Парциальное давление водяного пара [Па]
    Returns:
        Number: Атмосферное давление [Па]
    '''
    return e * (M + d) / d


@vec()
def rho_vT(rho: Number, vT: Number) -> Number:
    '''
    Расчет общего давления по плотности и виртуальной температуре
    Args:
        rho (Number): Плотность влажного воздуха [кг/м³]
        vT (Number): Виртуальная температура [K]
    Returns:
        Number: Атмосферное давление [Па]
    '''
    return R_air * vT * rho
