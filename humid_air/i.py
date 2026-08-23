'''Методы для работы с энтальпией'''

from .vec import vec, Number
from .consts import C_p_h2o, C_p_air, k1, k2


@vec()
def d_t(d: Number, t: Number) -> Number:
    '''
    Расчет энтальпии по влагосодержанию и температуре
    Args:
        d (Number): Влагосодержание [доля]
        t (Number): Температура воздуха [°C]
    Returns:
        Number: Энтальпия влажного воздуха [Дж/кг]
    '''
    return C_p_air * t + d * (k1 + (C_p_h2o - k2) * t)
