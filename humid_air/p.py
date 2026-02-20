'''Методы для работы с атмосферным давлением'''

from .vec import vec, ArrayLike
from .consts import CtoK, M, R_air


@vec(3)
def e_t_vT(e: ArrayLike, t: ArrayLike, vT: ArrayLike):
    '''
    Расчет общего давления по парциальному давлению, температуре и виртуальной температуре
    Args:
        e (ArrayLike): Парциальное давление водяного пара [Па]
        t (ArrayLike): Температура воздуха [°C]
        vT (ArrayLike): Виртуальная температура [K]
    Returns:
        ArrayLike: Атмосферное давление [Па]
    '''
    return (M * e * (t + CtoK)) / (vT - (t + CtoK))


@vec(2)
def d_e(d: ArrayLike, e: ArrayLike):
    '''
    Расчет общего давления по влагосодержанию и парциальному давлению
    Args:
        d (ArrayLike): Влагосодержание [доля]
        e (ArrayLike): Парциальное давление водяного пара [Па]
    Returns:
        ArrayLike: Атмосферное давление [Па]
    '''
    return e * (M + d) / d


@vec(2)
def rho_vT(rho: ArrayLike, vT: ArrayLike):
    '''
    Расчет общего давления по плотности и виртуальной температуре
    Args:
        rho (ArrayLike): Плотность влажного воздуха [кг/м³]
        vT (ArrayLike): Виртуальная температура [K]
    Returns:
        ArrayLike: Атмосферное давление [Па]
    '''
    return R_air * vT * rho
