'''Методы для работы с относительной влажностью'''

from .vec import vec, ArrayLike


@vec(2)
def E_sat_e_evap(E_sat: ArrayLike, e_evap: ArrayLike):
    '''
    Расчет относительной влажности по давлению насыщения и парциальному давлению
    Args:
        E_sat (ArrayLike): Давление насыщенного пара [Па]
        e_evap (ArrayLike): Парциальное давление водяного пара [Па]
    Returns:
        ArrayLike: Относительная влажность [доля]
    '''
    return e_evap / E_sat
