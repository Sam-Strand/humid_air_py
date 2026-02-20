'''Методы для работы с удельной теплотой парообразования'''

from .vec import vec, ArrayLike


@vec(2)
def d_L_sat(d: ArrayLike, L_sat: ArrayLike):
    '''
    Расчет общей теплоты парообразования
    Args:
        d (ArrayLike): Влагосодержание [доля]
        L_sat (ArrayLike): Удельная теплота образования пара [Дж/кг]
    Returns:
        ArrayLike: Общая теплота парообразования [Дж/кг]
    '''
    return L_sat * d
