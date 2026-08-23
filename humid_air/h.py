'''Методы для работы с относительной влажностью'''

from .vec import vec, Number


@vec()
def E_sat_e_evap(E_sat: Number, e_evap: Number) -> Number:
    '''
    Расчет относительной влажности по давлению насыщения и парциальному давлению
    Args:
        E_sat (Number): Давление насыщенного пара [Па]
        e_evap (Number): Парциальное давление водяного пара [Па]
    Returns:
        Number: Относительная влажность [доля]
    '''
    return e_evap / E_sat
