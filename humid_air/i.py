'''Методы для работы с энтальпией'''

from .vec import vec, ArrayLike
from .consts import C_p_h2o, C_p_air, k1, k2


@vec(2)
def d_t(d: ArrayLike, t: ArrayLike):
    '''
    Расчет энтальпии по влагосодержанию и температуре
    Args:
        d (ArrayLike): Влагосодержание [доля]
        t (ArrayLike): Температура воздуха [°C]
    Returns:
        ArrayLike: Энтальпия влажного воздуха [Дж/кг]
    '''
    return C_p_air * t + d * (k1 + (C_p_h2o - k2) * t)
