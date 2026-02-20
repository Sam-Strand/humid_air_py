'''Методы для работы с плотностью'''

from .vec import vec, ArrayLike
from .consts import R_air


@vec(2)
def p_vT(p: ArrayLike, vT: ArrayLike):
    '''
    Расчет плотности по давлению и виртуальной температуре
    Args:
        p (ArrayLike): Атмосферное давление [Па]
        vT (ArrayLike): Виртуальная температура [K]
    Returns:
        ArrayLike: Плотность влажного воздуха [кг/м³]
    '''
    return p / (R_air * vT)


@vec(1)
def V(V: ArrayLike):
    '''
    Расчет плотности по удельному объему
    Args:
        V (ArrayLike): Удельный объем [м³/кг]
    Returns:
        ArrayLike: Плотность влажного воздуха [кг/м³]
    '''
    return 1 / V
