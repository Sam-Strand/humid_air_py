'''
Модуль для расчета параметров влажного воздуха
E   [Па] - Давление насыщенного пара
e   [Па] - Парциальное давление водяного пара
t   [°C] - Температура воздуха
d   [доля] - Массовая доля водяного пара (влагосодержание)
rho [кг/м³] - Плотность влажного воздуха
p   [Па] - Атмосферное давление
h   [доля] - Относительная влажность
c   [Дж/(кг·K)] - Удельная теплоемкость влажного воздуха
i   [Дж/кг] - Энтальпия влажного воздуха
vT  [K] - Виртуальная температура
L   [Дж/кг] - Удельная теплота образования насыщенного пара
l   [Дж/кг] - Удельная теплота парообразования
v   [м/с] - скорость
V   [м³] - Объем
'''

from numba import vectorize
from math import exp, log
import numpy.typing as npt

# Универсальная газовая постоянная [Дж/(моль·K)]
R = 8.31446261815324
# Молярная масса воздуха [кг/моль]
M_air = 0.02898
# Молярная масса воды [кг/моль]
M_h2o = 0.01801528
# Отношение газовых постоянных [доля]
M = M_h2o / M_air
# Газовая постоянная для сухого воздуха [Дж/(кг·K)]
R_air = R / M_air
# Удельные теплоемкости сухого воздуха [Дж/(кг·K)]
C_p_air = 1005
# Удельные теплоемкости водяного пара [Дж/(кг·K)]
C_p_h2o = 1864
# Перевод [°C] в [K] и обратно
CtoK = 273.15
# Скрытая теплота парообразования [Дж]
k1 = 2501000
# Поправочный коэффициент температуры теплоты парообразования [Дж/K]
k2 = 2360


ArrayLike = npt.ArrayLike


def vec(n, nopython=True, cache=True):
    '''n аргументов float64, возвращает float64'''
    return vectorize([f"float64({', '.join(['float64'] * n)})"], nopython=nopython, cache=cache)


class Humid_air:
    '''
    Класс для расчета параметров влажного воздуха с поддержкой numpy
    '''
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
        E = Humid_air.E.t(t)
        e = Humid_air.e.E_h(E, h)
        return Humid_air.t.E(e)

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
        E = Humid_air.E.t(t)
        return Humid_air.d.e_p(E, p)

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
        E = Humid_air.E.t(t)
        e = Humid_air.e.E_h(E, h)
        vT = Humid_air.vT.e_p_t(e, p, t)
        return Humid_air.rho.p_vT(p, vT)

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
        E = Humid_air.E.t(t)
        e = Humid_air.e.E_h(E, h)
        return Humid_air.d.e_p(e, p)

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

    class L:
        '''Методы для работы с удельной теплотой образования пара'''

        @staticmethod
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

        @staticmethod
        @vec(2)
        def l_d(l: ArrayLike, d: ArrayLike):
            '''
            Расчет удельной теплоты на единицу влагосодержания

            Args:
                l (ArrayLike): Удельная теплота парообразования [Дж/кг]
                d (ArrayLike): Влагосодержание [доля]

            Returns:
                ArrayLike: Удельная теплота на единицу влагосодержания [Дж/кг]
            '''
            return l / d

    class l:
        '''Методы для работы с удельной теплотой парообразования'''

        @staticmethod
        @vec(2)
        def d_L(d: ArrayLike, L: ArrayLike):
            '''
            Расчет общей теплоты парообразования

            Args:
                d (ArrayLike): Влагосодержание [доля]
                L (ArrayLike): Удельная теплота образования пара [Дж/кг]

            Returns:
                ArrayLike: Общая теплота парообразования [Дж/кг]
            '''
            return L * d

    class E:
        '''Методы для работы с давлением насыщенного пара'''

        @staticmethod
        @vec(1)
        def t(t: ArrayLike):
            '''
            Расчет давления насыщенного пара по температуре (формула Магнуса)

            Args:
                t (ArrayLike): Температура воздуха [°C]

            Returns:
                ArrayLike: Давление насыщенного пара [Па]
            '''
            a = 611.2
            b = 17.62 if t > 0 else 22.46
            c = 243.12 if t > 0 else 272.62
            return a * exp(b * t / (c + t))

        @staticmethod
        @vec(2)
        def e_h(e: ArrayLike, h: ArrayLike):
            '''
            Расчет давления насыщенного пара по парциальному давлению и влажности

            Args:
                e (ArrayLike): Парциальное давление водяного пара [Па]
                h (ArrayLike): Относительная влажность [доля]

            Returns:
                ArrayLike: Давление насыщенного пара [Па]
            '''
            return e / h

    class e:
        '''Методы для работы с парциальным давлением водяного пара'''

        @staticmethod
        @vec(2)
        def E_h(E: ArrayLike, h: ArrayLike):
            '''
            Расчет парциального давления по давлению насыщения и влажности

            Args:
                E (ArrayLike): Давление насыщенного пара [Па]
                h (ArrayLike): Относительная влажность [доля]

            Returns:
                ArrayLike: Парциальное давление водяного пара [Па]
            '''
            return h * E

        @staticmethod
        @vec(2)
        def d_p(d: ArrayLike, p: ArrayLike):
            '''
            Расчет парциального давления по влагосодержанию и общему давлению

            Args:
                d (ArrayLike): Влагосодержание [доля]
                p (ArrayLike): Атмосферное давление [Па]

            Returns:
                ArrayLike: Парциальное давление водяного пара [Па]
            '''
            return (d * p) / (M + d)

        @staticmethod
        @vec(3)
        def p_t_vT(p: ArrayLike, t: ArrayLike, vT: ArrayLike):
            '''
            Расчет парциального давления по общему давлению, температуре и виртуальной температуре

            Args:
                p (ArrayLike): Атмосферное давление [Па]
                t (ArrayLike): Температура воздуха [°C]
                vT (ArrayLike): Виртуальная температура [K]

            Returns:
                ArrayLike: Парциальное давление водяного пара [Па]
            '''
            return (p * (vT - (t + CtoK))) / (M * (t + CtoK))

    class d:
        '''Методы для работы с влагосодержанием'''

        @staticmethod
        @vec(2)
        def e_p(e: ArrayLike, p: ArrayLike):
            '''
            Расчет влагосодержания по парциальному давлению и общему давлению

            Args:
                e (ArrayLike): Парциальное давление водяного пара [Па]
                p (ArrayLike): Атмосферное давление [Па]

            Returns:
                ArrayLike: Влагосодержание [доля]
            '''
            return M * e / (p - e)

        @staticmethod
        @vec(2)
        def L_l(L: ArrayLike, l: ArrayLike):
            '''
            Расчет влагосодержания по удельной и общей теплоте парообразования

            Args:
                L (ArrayLike): Удельная теплота образования пара [Дж/кг]
                l (ArrayLike): Удельная теплота парообразования [Дж/кг]

            Returns:
                ArrayLike: Влагосодержание [доля]
            '''
            return l / L

        @staticmethod
        @vec(1)
        def c(c: ArrayLike):
            '''
            Расчет влагосодержания по удельной теплоемкости

            Args:
                c (ArrayLike): Удельная теплоемкость влажного воздуха [Дж/(кг·K)]

            Returns:
                ArrayLike: Влагосодержание [доля]
            '''
            return (c - C_p_air) / (C_p_h2o - C_p_air)

        @staticmethod
        @vec(2)
        def i_t(i: ArrayLike, t: ArrayLike):
            '''
            Расчет влагосодержания по энтальпии и температуре

            Args:
                i (ArrayLike): Энтальпия влажного воздуха [Дж/кг]
                t (ArrayLike): Температура воздуха [°C]

            Returns:
                ArrayLike: Влагосодержание [доля]
            '''
            return (i - C_p_air * t) / (k1 + (C_p_h2o - k2) * t)

    class t:
        '''Методы для работы с температурой'''

        @staticmethod
        @vec(1)
        def L(L: ArrayLike):
            '''
            Расчет температуры по удельной теплоте образования пара

            Args:
                L (ArrayLike): Удельная теплота образования пара [Дж/кг]

            Returns:
                ArrayLike: Температура [°C]
            '''
            return (k1 + L) / k2

        @staticmethod
        @vec(3)
        def e_p_vT(e: ArrayLike, p: ArrayLike, vT: ArrayLike):
            '''
            Расчет температуры по парциальному давлению, общему давлению и виртуальной температуре

            Args:
                e (ArrayLike): Парциальное давление водяного пара [Па]
                p (ArrayLike): Атмосферное давление [Па]
                vT (ArrayLike): Виртуальная температура [K]

            Returns:
                ArrayLike: Температура воздуха [°C]
            '''
            return (vT / (1 + M * e / p)) - CtoK

        @staticmethod
        @vec(1)
        def E(E: ArrayLike):
            '''
            Расчет температуры по давлению насыщенного пара (обратная формула Магнуса)

            Args:
                E (ArrayLike): Давление насыщенного пара [Па]

            Returns:
                ArrayLike: Температура воздуха [°C]
            '''
            a = 611.2
            b = 17.62
            c = 243.12
            numerator = log(E / a)

            if numerator < 0:
                b = 22.46
                c = 272.62
                numerator = log(E / a)

            return c * numerator / (b - numerator)

        @staticmethod
        @vec(2)
        def i_d(i: ArrayLike, d: ArrayLike):
            '''
            Расчет температуры по энтальпии и влагосодержанию

            Args:
                i (ArrayLike): Энтальпия влажного воздуха [Дж/кг]
                d (ArrayLike): Влагосодержание [доля]

            Returns:
                ArrayLike: Температура воздуха [°C]
            '''
            return (i - k1 * d) / (C_p_air + d * (C_p_h2o - k2))

    class rho:
        '''Методы для работы с плотностью'''

        @staticmethod
        @vec(2)
        def p_vT(p: ArrayLike, vT: ArrayLike):
            '''
            Расчет плотности по давлению и виртуальной температуре

            Args:
                p (ArrayLike): Атмосферное давление [Па]
                vT (ArrayLike): Виртуальная температура [K]

            Returns:
                ArrayLike: Плотность влажного воздуха [кг/м³]
            '''
            return p / (R_air * vT)

        @staticmethod
        @vec(1)
        def V(V: ArrayLike):
            '''
            Расчет плотности по удельному объему

            Args:
                V (ArrayLike): Удельный объем [м³/кг]

            Returns:
                ArrayLike: Плотность влажного воздуха [кг/м³]
            '''
            return 1 / V

    class h:
        '''Методы для работы с относительной влажностью'''

        @staticmethod
        @vec(2)
        def E_e(E: ArrayLike, e: ArrayLike):
            '''
            Расчет относительной влажности по давлению насыщения и парциальному давлению

            Args:
                E (ArrayLike): Давление насыщенного пара [Па]
                e (ArrayLike): Парциальное давление водяного пара [Па]

            Returns:
                ArrayLike: Относительная влажность [доля]
            '''
            return e / E

    class i:
        '''Методы для работы с энтальпией'''

        @staticmethod
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

    class vT:
        '''Методы для работы с виртуальной температурой'''

        @staticmethod
        @vec(3)
        def e_p_t(e: ArrayLike, p: ArrayLike, t: ArrayLike):
            '''
            Расчет виртуальной температуры по парциальному давлению, общему давлению и температуре

            Args:
                e (ArrayLike): Парциальное давление водяного пара [Па]
                p (ArrayLike): Атмосферное давление [Па]
                t (ArrayLike): Температура воздуха [°C]

            Returns:
                ArrayLike: Виртуальная температура [K]
            '''
            return (t + CtoK) / (1 - (1 - M) * e / p)

        @staticmethod
        @vec(2)
        def p_rho(p: ArrayLike, rho: ArrayLike):
            '''
            Расчет виртуальной температуры по давлению и плотности

            Args:
                p (ArrayLike): Атмосферное давление [Па]
                rho (ArrayLike): Плотность влажного воздуха [кг/м³]

            Returns:
                ArrayLike: Виртуальная температура [K]
            '''
            return p / (R_air * rho)

    class p:
        '''Методы для работы с атмосферным давлением'''

        @staticmethod
        @vec(3)
        def e_t_vT(e: ArrayLike, t: ArrayLike, vT: ArrayLike):
            '''
            Расчет общего давления по парциальному давлению, температуре и виртуальной температуре

            Args:
                e (ArrayLike): Парциальное давление водяного пара [Па]
                t (ArrayLike): Температура воздуха [°C]
                vT (ArrayLike): Виртуальная температура [K]

            Returns:
                ArrayLike: Атмосферное давление [Па]
            '''
            return (M * e * (t + CtoK)) / (vT - (t + CtoK))

        @staticmethod
        @vec(2)
        def d_e(d: ArrayLike, e: ArrayLike):
            '''
            Расчет общего давления по влагосодержанию и парциальному давлению

            Args:
                d (ArrayLike): Влагосодержание [доля]
                e (ArrayLike): Парциальное давление водяного пара [Па]

            Returns:
                ArrayLike: Атмосферное давление [Па]
            '''
            return e * (M + d) / d

        @staticmethod
        @vec(2)
        def rho_vT(rho: ArrayLike, vT: ArrayLike):
            '''
            Расчет общего давления по плотности и виртуальной температуре

            Args:
                rho (ArrayLike): Плотность влажного воздуха [кг/м³]
                vT (ArrayLike): Виртуальная температура [K]

            Returns:
                ArrayLike: Атмосферное давление [Па]
            '''
            return R_air * vT * rho

    class c:
        '''Методы для работы с удельной теплоемкостью'''

        @staticmethod
        @vec(1)
        def d(d: ArrayLike):
            '''
            Расчет удельной теплоемкости по влагосодержанию

            Args:
                d (ArrayLike): Влагосодержание [доля]

            Returns:
                ArrayLike: Удельная теплоемкость влажного воздуха [Дж/(кг·K)]
            '''
            return C_p_air * (1 - d) + C_p_h2o * d

    class V:
        '''Методы для работы с удельным объемом'''

        @staticmethod
        @vec(1)
        def rho(rho: ArrayLike):
            '''
            Расчет удельного объема по плотности

            Args:
                rho (ArrayLike): Плотность влажного воздуха [кг/м³]

            Returns:
                ArrayLike: Удельный объем [м³/кг]
            '''
            return 1 / rho


if __name__ == '__main__':
    # Пример расчетов
    print(Humid_air.density(10, 101325, 0.5))  # 1.244425737661975
    print(Humid_air.density([10], (101325), [0.5]))  # [1.24442581]
    # [1.2467098  1.24613878 1.24556777 1.24499675 1.24442574 1.24385472 1.24328371 1.24271269 1.24214168 1.24157066]
    print(Humid_air.density(10, 101325,
                            [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]))
