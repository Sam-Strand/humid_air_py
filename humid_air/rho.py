'''Методы для работы с плотностью'''

from .vec import vec, Number
from .consts import R_air


@vec()
def p_vT(p: Number, vT: Number) -> Number:
    '''
    Расчет плотности по давлению и виртуальной температуре
    Args:
        p (Number): Атмосферное давление [Па]
        vT (Number): Виртуальная температура [K]
    Returns:
        Number: Плотность влажного воздуха [кг/м³]
    '''
    return p / (R_air * vT)


@vec()
def V(V: Number) -> Number:
    '''
    Расчет плотности по удельному объему
    Args:
        V (Number): Удельный объем [м³/кг]
    Returns:
        Number: Плотность влажного воздуха [кг/м³]
    '''
    return 1 / V
