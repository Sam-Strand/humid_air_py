'''Методы для работы с удельной теплоемкостью'''

from .vec import vec, ArrayLike
from .consts import C_p_air, C_p_h2o


@vec(1)
def d(d: ArrayLike):
    '''
    Расчет удельной теплоемкости по влагосодержанию
    Args:
        d (ArrayLike): Влагосодержание [доля]
    Returns:
        ArrayLike: Удельная теплоемкость влажного воздуха [Дж/(кг·K)]
    '''
    return C_p_air * (1 - d) + C_p_h2o * d
