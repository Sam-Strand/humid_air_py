'''Методы для работы с удельной теплотой парообразования'''

from .vec import vec, Number


@vec()
def d_L_sat(d: Number, L_sat: Number) -> Number:
    '''
    Расчет общей теплоты парообразования
    Args:
        d (Number): Влагосодержание [доля]
        L_sat (Number): Удельная теплота образования пара [Дж/кг]
    Returns:
        Number: Общая теплота парообразования [Дж/кг]
    '''
    return L_sat * d
