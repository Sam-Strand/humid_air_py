'''Методы для работы с удельной теплотой образования пара'''

from .vec import vec, Number
from .consts import k1, k2

@vec()
def t(t: Number) -> Number:
    '''
    Расчет удельной теплоты образования пара по температуре
    Args:
        t (Number): Температура [°C]
    Returns:
        Number: Удельная теплота образования пара [Дж/кг]
    '''
    return k1 - k2 * t


@vec()
def l_evap_d(l_evap: Number, d: Number) -> Number:
    '''
    Расчет удельной теплоты на единицу влагосодержания
    Args:
        l_evap (Number): Удельная теплота парообразования [Дж/кг]
        d (Number): Влагосодержание [доля]
    Returns:
        Number: Удельная теплота на единицу влагосодержания [Дж/кг]
    '''
    return l_evap / d
