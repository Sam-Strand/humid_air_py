'''Методы для работы с удельным объемом'''

from .vec import vec, Number


@vec()
def rho(rho: Number) -> Number:
    '''
    Расчет удельного объема по плотности
    Args:
        rho (Number): Плотность влажного воздуха [кг/м³]
    Returns:
        Number: Удельный объем [м³/кг]
    '''
    return 1 / rho
