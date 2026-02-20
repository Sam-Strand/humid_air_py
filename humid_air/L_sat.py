'''Методы для работы с удельной теплотой образования пара'''

from .vec import vec, ArrayLike
from .consts import k1, k2

@vec(1)
def t(t: ArrayLike):
    '''
    Расчет удельной теплоты образования пара по температуре
    Args:
        t (ArrayLike): Температура [°C]
    Returns:
        ArrayLike: Удельная теплота образования пара [Дж/кг]
    '''
    return k1 - k2 * t


@vec(2)
def l_evap_d(l_evap: ArrayLike, d: ArrayLike):
    '''
    Расчет удельной теплоты на единицу влагосодержания
    Args:
        l_evap (ArrayLike): Удельная теплота парообразования [Дж/кг]
        d (ArrayLike): Влагосодержание [доля]
    Returns:
        ArrayLike: Удельная теплота на единицу влагосодержания [Дж/кг]
    '''
    return l_evap / d