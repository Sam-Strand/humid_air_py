from math import exp, log
from .vec import vec, ArrayLike
from .consts import C_p_water

from .E_sat import t as E_sat_t
from .t import i_d as t_i_d, E_sat as t_E_sat
from .e_evap import E_sat_h as e_evap_E_sat_h, d_p as e_evap_d_p
from .d import e_evap_p as d_e_evap_p
from .i import d_t as i_d_t
from .h import E_sat_e_evap as h_E_sat_e_evap

from .vT import e_p_t as vT_e_p_t
from .rho import p_vT as rho_p_vT

class Humid_air:
    '''
    Класс для расчета параметров влажного воздуха с поддержкой numpy
    '''
    @staticmethod
    @vec(6)
    def g_water(t: ArrayLike, p: ArrayLike, h: ArrayLike, t_water: ArrayLike, minT: ArrayLike, maxH: ArrayLike):
        '''
        Оптимизация удельного расхода воды [кг].
        Возвращает удельный расход.
        '''
        if h >= 1.0:
            return 0.0
        E_sat = E_sat_t(t)
        
        e = e_evap_E_sat_h(E_sat, h)
        d = d_e_evap_p(e, p)
        d_max = d_e_evap_p(E_sat, p)
        i = i_d_t(d, t)
        i_w = C_p_water * t_water
        max_specific_water = max(0, (d_max - d) * 1.2)

        low_g = 0.0
        high_g = max_specific_water
        g = 0.0

        for _ in range(50):
            mid_g = (low_g + high_g) / 2
            
            d2 = d + mid_g
            i2 = i + mid_g * i_w
            t_new = t_i_d(i2, d2)
            E2 = E_sat_t(t_new)
            e2 = e_evap_d_p(d2, p)
            h_new = h_E_sat_e_evap(E2, e2)
            
            conditions_met = (t_new >= minT) and (h_new <= maxH)
            
            if conditions_met:
                g = mid_g
                low_g = mid_g
            else:
                high_g = mid_g

            if abs(high_g - low_g) < 1e-6:
                break

        return g

    @staticmethod
    @vec(4)
    def heating_target_temperature(p: ArrayLike, t1: ArrayLike, h1: ArrayLike, h2: ArrayLike) -> ArrayLike:
        '''
        Расчет температуры нагрева для достижения целевой влажности при постоянном влагосодержании

        Процесс нагрева происходит при d=const (изохорный нагрев). При переходе через 0°C 
        используются разные коэффициенты формулы Магнуса для льда и воды.

        Args:
            p (ArrayLike): Атмосферное давление [Па]
            t1 (ArrayLike): Начальная температура воздуха [°C]
            h1 (ArrayLike): Начальная относительная влажность [доля]
            h2 (ArrayLike): Целевая относительная влажность [доля]

        Returns:
            ArrayLike: Температура [°C], при которой влажность достигает h2

        Алгоритм:
        1. Если начальная температура >= 0°C:
        - Используются коэффициенты формулы Магнуса для воды
        - Прямой расчет конечной температуры

        2. Если начальная температура < 0°C:
        - Рассчитывается влажность при 0°C
        - Если h(0°C) <= h2: цель достигается ДО 0°C (используются коэффициенты для льда)
        - Если h(0°C) > h2: цель достигается ПОСЛЕ 0°C (используются коэффициенты для воды)
        '''
        a = 611.2
        b_water, c_water = 17.62, 243.12
        if t1 >= 0:
            E1 = a * exp(b_water * t1 / (c_water + t1))
            d = 0.622 * (h1 * E1) / (p - h1 * E1)
            E2 = (d * p) / (0.622 * h2 + d * h2)
            K = log(E2 / a)
            t2 = c_water * K / (b_water - K)
            return t2
        else:
            b_ice, c_ice = 22.46, 272.62
            E1 = a * exp(b_ice * t1 / (c_ice + t1))
            d = 0.622 * (h1 * E1) / (p - h1 * E1)
            E_at_0_water = a * exp(b_water * 0 / (c_water + 0))
            h_at_0 = (d * p) / (0.622 * E_at_0_water + d * E_at_0_water)
            E2 = (d * p) / (0.622 * h2 + d * h2)
            K = log(E2 / a)
            if h_at_0 <= h2:
                return c_ice * K / (b_ice - K)
            else:
                return c_water * K / (b_water - K)

    @staticmethod
    def dew_point_temperature(t: ArrayLike, h: ArrayLike):
        '''
        Расчет температуры точки росы по температуре и относительной влажности

        Args:
            t (ArrayLike): Температура воздуха [°C]
            h (ArrayLike): Относительная влажность [доля]

        Returns:
            ArrayLike: Температура точки росы [°C]
        '''
        E_sat = E_sat_t(t)
        e = e_evap_E_sat_h(E_sat, h)
        return t_E_sat(e)

    @staticmethod
    def maximum_moisture_content(p: ArrayLike, t: ArrayLike):
        '''
        Расчет предельного влагосодержания (при 100% влажности)

        Args:
            p (ArrayLike): Атмосферное давление [Па]
            t (ArrayLike): Температура воздуха [°C]

        Returns:
            ArrayLike: Предельное влагосодержание [доля]
        '''
        E_sat = E_sat_t(t)
        return d_e_evap_p(E_sat, p)

    @staticmethod
    def density(t: ArrayLike, p: ArrayLike, h: ArrayLike):
        '''
        Расчет плотности влажного воздуха

        Args:
            t (ArrayLike): Температура воздуха [°C]
            p (ArrayLike): Атмосферное давление [Па]
            h (ArrayLike): Относительная влажность [доля]

        Returns:
            ArrayLike: Плотность влажного воздуха [кг/м³]
        '''
        E_sat = E_sat_t(t)
        e = e_evap_E_sat_h(E_sat, h)
        vT = vT_e_p_t(e, p, t)
        return rho_p_vT(p, vT)

    @staticmethod
    def moisture_content(t: ArrayLike, p: ArrayLike, h: ArrayLike):
        '''
        Расчет влагосодержания по температуре, давлению и влажности

        Args:
            t (ArrayLike): Температура воздуха [°C]
            p (ArrayLike): Атмосферное давление [Па]
            h (ArrayLike): Относительная влажность [доля]

        Returns:
            ArrayLike: Влагосодержание [доля]
        '''
        E_sat = E_sat_t(t)
        e = e_evap_E_sat_h(E_sat, h)
        return d_e_evap_p(e, p)

    @staticmethod
    @vec(2)
    def dynamic_pressure(rho: ArrayLike, v: ArrayLike):
        '''
        Расчет динамического давления

        Args:
            rho (ArrayLike): Плотность воздуха [кг/м³]
            v (ArrayLike): Скорость воздуха [м/с]

        Returns:
            ArrayLike: Динамическое давление [Па]
        '''
        return rho * v ** 2 / 2

    class L_sat:
        from .L_sat import t, l_evap_d

    class l_evap:
        from .l_evap import d_L_sat

    class E_sat:
        from .E_sat import e_evap_h, t

    class e_evap:
        from .e_evap import E_sat_h, d_p, p_t_vT

    class d:
        from .d import e_evap_p, L_l, c, i_t
    class t:
        from .t import L, e_evap_p_vT, E_sat, i_d

    class rho:
        from .rho import p_vT, V

    class h:
        from .h import E_sat_e_evap

    class i:
        from .i import d_t

    class vT:
        from .vT import e_p_t, p_rho

    class p:
        from .p import e_t_vT, d_e, rho_vT

    class c:
        from .c import d

    class V:
        from .V import rho
