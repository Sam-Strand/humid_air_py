'''Методы для работы с удельной теплоемкостью'''

from .vec import vec, Number
from .consts import C_p_air, C_p_h2o


@vec()
def d(d: Number) -> Number:
    '''
    Расчет удельной теплоемкости по влагосодержанию
    Args:
        d (Number): Влагосодержание [доля]
    Returns:
        Number: Удельная теплоемкость влажного воздуха [Дж/(кг·K)]
    '''
    return C_p_air * (1 - d) + C_p_h2o * d
