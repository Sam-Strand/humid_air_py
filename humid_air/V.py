'''Методы для работы с удельным объемом'''

from .vec import vec, ArrayLike


@vec(1)
def rho(rho: ArrayLike):
    '''
    Расчет удельного объема по плотности
    Args:
        rho (ArrayLike): Плотность влажного воздуха [кг/м³]
    Returns:
        ArrayLike: Удельный объем [м³/кг]
    '''
    return 1 / rho
