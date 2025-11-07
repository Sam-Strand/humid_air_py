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

import numpy as np
from numba import vectorize, float64

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


class Humid_air:
    '''
    Класс для расчета параметров влажного воздуха с поддержкой numpy
    '''
    @vectorize([float64(float64, float64, float64, float64)], nopython=True, cache=True)
    def heating_target_temperature(p, t1, h1, h2):
        '''
        Расчет температуры нагрева для достижения целевой влажности при постоянном влагосодержании
        
        Процесс нагрева происходит при d=const (изохорный нагрев). При переходе через 0°C 
        используются разные коэффициенты формулы Магнуса для льда и воды.
        
        Args:
            p (float): Атмосферное давление [Па]
            t1 (float): Начальная температура воздуха [°C]
            h1 (float): Начальная относительная влажность [доля]
            h2 (float): Целевая относительная влажность [доля]
            
        Returns:
            float: Температура [°C], при которой влажность достигает h2
            
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
        b_ice, c_ice = 22.46, 272.62

        if t1 >= 0:
            E1 = a * np.exp(b_water * t1 / (c_water + t1))
            d = 0.622 * (h1 * E1) / (p - h1 * E1)
            E2 = (d * p) / (0.622 * h2 + d * h2)
            K = np.log(E2 / a)
            t2 = c_water * K / (b_water - K)
            return t2

        else:
            E1_ice = a * np.exp(b_ice * t1 / (c_ice + t1))
            d = 0.622 * (h1 * E1_ice) / (p - h1 * E1_ice)

            E_at_0_water = a * np.exp(b_water * 0 / (c_water + 0))
            h_at_0 = (d * p) / (0.622 * E_at_0_water + d * E_at_0_water)

            if h_at_0 <= h2:
                E2 = (d * p) / (0.622 * h2 + d * h2)
                K = np.log(E2 / a)
                t2 = c_ice * K / (b_ice - K)
                return t2

            else:
                E2 = (d * p) / (0.622 * h2 + d * h2)
                K = np.log(E2 / a)
                t2 = c_water * K / (b_water - K)
                return t2

    @staticmethod
    @vectorize([float64(float64)], nopython=True, cache=True)
    def magnus_pressure(t):
        '''
        Расчет давления насыщенного пара по температуре (формула Магнуса)

        Args:
            t (float или np.ndarray): Температура воздуха [°C]

        Returns:
            float или np.ndarray: Давление насыщенного пара [Па]
        '''
        a = 611.2
        b = 17.62 if t > 0 else 22.46
        c = 243.12 if t > 0 else 272.62
        return a * np.exp(b * t / (c + t))

    @staticmethod
    @vectorize([float64(float64)], nopython=True, cache=True)
    def inverse_magnus(E):
        '''
        Расчет температуры по давлению насыщенного пара (обратная формула Магнуса)

        Args:
            E (float или np.ndarray): Давление насыщенного пара [Па]

        Returns:
            float или np.ndarray: Температура воздуха [°C]
        '''
        a = 611.2
        # Пробуем сначала для положительных температур
        b = 17.62
        c = 243.12
        numerator = np.log(E / a)

        # Если не получилось (отрицательная температура), пробуем для отрицательных
        if numerator < 0:
            b = 22.46
            c = 272.62
            numerator = np.log(E / a)

        return c * numerator / (b - numerator)

    @staticmethod
    def dew_point_temperature(t, h):
        '''
        Расчет температуры точки росы по температуре и относительной влажности

        Args:
            t (float или np.ndarray): Температура воздуха [°C]
            h (float или np.ndarray): Относительная влажность [доля]

        Returns:
            float или np.ndarray: Температура точки росы [°C]
        '''
        E = Humid_air.E.t(t)
        e = Humid_air.e.E_h(E, h)
        return Humid_air.t.E(e)

    @staticmethod
    def maximum_moisture_content(p, t):
        '''
        Расчет предельного влагосодержания (при 100% влажности)

        Args:
            p (float или np.ndarray): Атмосферное давление [Па]
            t (float или np.ndarray): Температура воздуха [°C]

        Returns:
            float или np.ndarray: Предельное влагосодержание [доля]
        '''
        E = Humid_air.E.t(t)
        return Humid_air.d.e_p(E, p)

    @staticmethod
    def density(t, p, h):
        '''
        Расчет плотности влажного воздуха

        Args:
            t (float или np.ndarray): Температура воздуха [°C]
            p (float или np.ndarray): Атмосферное давление [Па]
            h (float или np.ndarray): Относительная влажность [доля]

        Returns:
            float или np.ndarray: Плотность влажного воздуха [кг/м³]
        '''
        E = Humid_air.E.t(t)
        e = Humid_air.e.E_h(E, h)
        vT = Humid_air.vT.e_p_t(e, p, t)
        return Humid_air.rho.p_vT(p, vT)

    @staticmethod
    def moisture_content(t, p, h):
        '''
        Расчет влагосодержания по температуре, давлению и влажности

        Args:
            t (float или np.ndarray): Температура воздуха [°C]
            p (float или np.ndarray): Атмосферное давление [Па]
            h (float или np.ndarray): Относительная влажность [доля]

        Returns:
            float или np.ndarray: Влагосодержание [доля]
        '''
        E = Humid_air.E.t(t)
        e = Humid_air.e.E_h(E, h)
        return Humid_air.d.e_p(e, p)

    @staticmethod
    def dynamic_pressure(rho, v):
        '''
        Расчет динамического давления

        Args:
            rho (float или np.ndarray): Плотность воздуха [кг/м³]
            v (float или np.ndarray): Скорость воздуха [м/с]

        Returns:
            float или np.ndarray: Динамическое давление [Па]
        '''
        rho = np.asarray(rho, dtype=np.float64)
        v = np.asarray(v, dtype=np.float64)
        return rho * v ** 2 / 2

    class L:
        '''Методы для работы с удельной теплотой образования пара'''

        @staticmethod
        def t(t):
            '''
            Расчет удельной теплоты образования пара по температуре

            Args:
                t (float или np.ndarray): Температура [°C]

            Returns:
                float или np.ndarray: Удельная теплота образования пара [Дж/кг]
            '''
            t = np.asarray(t, dtype=np.float64)
            return k1 - k2 * t

        @staticmethod
        def l_d(l, d):
            '''
            Расчет удельной теплоты на единицу влагосодержания

            Args:
                l (float или np.ndarray): Удельная теплота парообразования [Дж/кг]
                d (float или np.ndarray): Влагосодержание [доля]

            Returns:
                float или np.ndarray: Удельная теплота на единицу влагосодержания [Дж/кг]
            '''
            l = np.asarray(l, dtype=np.float64)
            d = np.asarray(d, dtype=np.float64)
            return l / d

    class l:
        '''Методы для работы с удельной теплотой парообразования'''

        @staticmethod
        def d_L(d, L):
            '''
            Расчет общей теплоты парообразования

            Args:
                d (float или np.ndarray): Влагосодержание [доля]
                L (float или np.ndarray): Удельная теплота образования пара [Дж/кг]

            Returns:
                float или np.ndarray: Общая теплота парообразования [Дж/кг]
            '''
            L = np.asarray(L, dtype=np.float64)
            d = np.asarray(d, dtype=np.float64)
            return L * d

    class E:
        '''Методы для работы с давлением насыщенного пара'''

        @staticmethod
        def t(t):
            '''
            Расчет давления насыщенного пара по температуре (формула Магнуса)

            Args:
                t (float или np.ndarray): Температура воздуха [°C]

            Returns:
                float или np.ndarray: Давление насыщенного пара [Па]
            '''
            return Humid_air.magnus_pressure(t)

        @staticmethod
        def e_h(e, h):
            '''
            Расчет давления насыщенного пара по парциальному давлению и влажности

            Args:
                e (float или np.ndarray): Парциальное давление водяного пара [Па]
                h (float или np.ndarray): Относительная влажность [доля]

            Returns:
                float или np.ndarray: Давление насыщенного пара [Па]
            '''
            e = np.asarray(e, dtype=np.float64)
            h = np.asarray(h, dtype=np.float64)
            return e / h

    class e:
        '''Методы для работы с парциальным давлением водяного пара'''

        @staticmethod
        def E_h(E, h):
            '''
            Расчет парциального давления по давлению насыщения и влажности

            Args:
                E (float или np.ndarray): Давление насыщенного пара [Па]
                h (float или np.ndarray): Относительная влажность [доля]

            Returns:
                float или np.ndarray: Парциальное давление водяного пара [Па]
            '''
            E = np.asarray(E, dtype=np.float64)
            h = np.asarray(h, dtype=np.float64)
            return h * E

        @staticmethod
        def d_p(d, p):
            '''
            Расчет парциального давления по влагосодержанию и общему давлению

            Args:
                d (float или np.ndarray): Влагосодержание [доля]
                p (float или np.ndarray): Атмосферное давление [Па]

            Returns:
                float или np.ndarray: Парциальное давление водяного пара [Па]
            '''
            d = np.asarray(d, dtype=np.float64)
            p = np.asarray(p, dtype=np.float64)
            return (d * p) / (M + d)

        @staticmethod
        def p_t_vT(p, t, vT):
            '''
            Расчет парциального давления по общему давлению, температуре и виртуальной температуре

            Args:
                p (float или np.ndarray): Атмосферное давление [Па]
                t (float или np.ndarray): Температура воздуха [°C]
                vT (float или np.ndarray): Виртуальная температура [K]

            Returns:
                float или np.ndarray: Парциальное давление водяного пара [Па]
            '''
            p = np.asarray(p, dtype=np.float64)
            t = np.asarray(t, dtype=np.float64)
            vT = np.asarray(vT, dtype=np.float64)
            return (p * (vT - (t + CtoK))) / (M * (t + CtoK))

    class d:
        '''Методы для работы с влагосодержанием'''

        @staticmethod
        def e_p(e, p):
            '''
            Расчет влагосодержания по парциальному давлению и общему давлению

            Args:
                e (float или np.ndarray): Парциальное давление водяного пара [Па]
                p (float или np.ndarray): Атмосферное давление [Па]

            Returns:
                float или np.ndarray: Влагосодержание [доля]
            '''
            e = np.asarray(e, dtype=np.float64)
            p = np.asarray(p, dtype=np.float64)
            return M * e / (p - e)

        @staticmethod
        def L_l(L, l):
            '''
            Расчет влагосодержания по удельной и общей теплоте парообразования

            Args:
                L (float или np.ndarray): Удельная теплота образования пара [Дж/кг]
                l (float или np.ndarray): Удельная теплота парообразования [Дж/кг]

            Returns:
                float или np.ndarray: Влагосодержание [доля]
            '''
            l = np.asarray(l, dtype=np.float64)
            L = np.asarray(L, dtype=np.float64)
            return l / L

        @staticmethod
        def c(c):
            '''
            Расчет влагосодержания по удельной теплоемкости

            Args:
                c (float или np.ndarray): Удельная теплоемкость влажного воздуха [Дж/(кг·K)]

            Returns:
                float или np.ndarray: Влагосодержание [доля]
            '''
            c = np.asarray(c, dtype=np.float64)
            return (c - C_p_air) / (C_p_h2o - C_p_air)

        @staticmethod
        def i_t(i, t):
            '''
            Расчет влагосодержания по энтальпии и температуре

            Args:
                i (float или np.ndarray): Энтальпия влажного воздуха [Дж/кг]
                t (float или np.ndarray): Температура воздуха [°C]

            Returns:
                float или np.ndarray: Влагосодержание [доля]
            '''
            i = np.asarray(i, dtype=np.float64)
            t = np.asarray(t, dtype=np.float64)
            return (i - C_p_air * t) / (k1 + (C_p_h2o - k2) * t)

    class t:
        '''Методы для работы с температурой'''

        @staticmethod
        def L(L):
            '''
            Расчет температуры по удельной теплоте образования пара

            Args:
                L (float или np.ndarray): Удельная теплота образования пара [Дж/кг]

            Returns:
                float или np.ndarray: Температура [°C]
            '''
            L = np.asarray(L, dtype=np.float64)
            return (k1 + L) / k2

        @staticmethod
        def e_p_vT(e, p, vT):
            '''
            Расчет температуры по парциальному давлению, общему давлению и виртуальной температуре

            Args:
                e (float или np.ndarray): Парциальное давление водяного пара [Па]
                p (float или np.ndarray): Атмосферное давление [Па]
                vT (float или np.ndarray): Виртуальная температура [K]

            Returns:
                float или np.ndarray: Температура воздуха [°C]
            '''
            e = np.asarray(e, dtype=np.float64)
            p = np.asarray(p, dtype=np.float64)
            vT = np.asarray(vT, dtype=np.float64)
            return (vT / (1 + M * e / p)) - CtoK

        @staticmethod
        def E(E):
            '''
            Расчет температуры по давлению насыщенного пара (обратная формула Магнуса)

            Args:
                E (float или np.ndarray): Давление насыщенного пара [Па]

            Returns:
                float или np.ndarray: Температура воздуха [°C]
            '''
            E = np.asarray(E, dtype=np.float64)
            return Humid_air.inverse_magnus(E)

        @staticmethod
        def i_d(i, d):
            '''
            Расчет температуры по энтальпии и влагосодержанию

            Args:
                i (float или np.ndarray): Энтальпия влажного воздуха [Дж/кг]
                d (float или np.ndarray): Влагосодержание [доля]

            Returns:
                float или np.ndarray: Температура воздуха [°C]
            '''
            i = np.asarray(i, dtype=np.float64)
            d = np.asarray(d, dtype=np.float64)
            return (i - k1 * d) / (C_p_air + d * (C_p_h2o - k2))

    class rho:
        '''Методы для работы с плотностью'''

        @staticmethod
        def p_vT(p, vT):
            '''
            Расчет плотности по давлению и виртуальной температуре

            Args:
                p (float или np.ndarray): Атмосферное давление [Па]
                vT (float или np.ndarray): Виртуальная температура [K]

            Returns:
                float или np.ndarray: Плотность влажного воздуха [кг/м³]
            '''
            p = np.asarray(p, dtype=np.float64)
            vT = np.asarray(vT, dtype=np.float64)
            return p / (R_air * vT)

        @staticmethod
        def V(V):
            '''
            Расчет плотности по удельному объему

            Args:
                V (float или np.ndarray): Удельный объем [м³/кг]

            Returns:
                float или np.ndarray: Плотность влажного воздуха [кг/м³]
            '''
            V = np.asarray(V, dtype=np.float64)
            return 1 / V

    class h:
        '''Методы для работы с относительной влажностью'''

        @staticmethod
        def E_e(E, e):
            '''
            Расчет относительной влажности по давлению насыщения и парциальному давлению

            Args:
                E (float или np.ndarray): Давление насыщенного пара [Па]
                e (float или np.ndarray): Парциальное давление водяного пара [Па]

            Returns:
                float или np.ndarray: Относительная влажность [доля]
            '''
            E = np.asarray(E, dtype=np.float64)
            e = np.asarray(e, dtype=np.float64)
            return e / E

    class i:
        '''Методы для работы с энтальпией'''

        @staticmethod
        def d_t(d, t):
            '''
            Расчет энтальпии по влагосодержанию и температуре

            Args:
                d (float или np.ndarray): Влагосодержание [доля]
                t (float или np.ndarray): Температура воздуха [°C]

            Returns:
                float или np.ndarray: Энтальпия влажного воздуха [Дж/кг]
            '''
            t = np.asarray(t, dtype=np.float64)
            d = np.asarray(d, dtype=np.float64)
            return C_p_air * t + d * (k1 + (C_p_h2o - k2) * t)

    class vT:
        '''Методы для работы с виртуальной температурой'''

        @staticmethod
        def e_p_t(e, p, t):
            '''
            Расчет виртуальной температуры по парциальному давлению, общему давлению и температуре

            Args:
                e (float или np.ndarray): Парциальное давление водяного пара [Па]
                p (float или np.ndarray): Атмосферное давление [Па]
                t (float или np.ndarray): Температура воздуха [°C]

            Returns:
                float или np.ndarray: Виртуальная температура [K]
            '''
            e = np.asarray(e, dtype=np.float64)
            p = np.asarray(p, dtype=np.float64)
            t = np.asarray(t, dtype=np.float64)
            return (t + CtoK) / (1 - (1 - M) * e / p)

        @staticmethod
        def p_rho(p, rho):
            '''
            Расчет виртуальной температуры по давлению и плотности

            Args:
                p (float или np.ndarray): Атмосферное давление [Па]
                rho (float или np.ndarray): Плотность влажного воздуха [кг/м³]

            Returns:
                float или np.ndarray: Виртуальная температура [K]
            '''
            p = np.asarray(p, dtype=np.float64)
            rho = np.asarray(rho, dtype=np.float64)
            return p / (R_air * rho)

    class p:
        '''Методы для работы с атмосферным давлением'''

        @staticmethod
        def e_t_vT(e, t, vT):
            '''
            Расчет общего давления по парциальному давлению, температуре и виртуальной температуре

            Args:
                e (float или np.ndarray): Парциальное давление водяного пара [Па]
                t (float или np.ndarray): Температура воздуха [°C]
                vT (float или np.ndarray): Виртуальная температура [K]

            Returns:
                float или np.ndarray: Атмосферное давление [Па]
            '''
            e = np.asarray(e, dtype=np.float64)
            t = np.asarray(t, dtype=np.float64)
            vT = np.asarray(vT, dtype=np.float64)
            return (M * e * (t + CtoK)) / (vT - (t + CtoK))

        @staticmethod
        def d_e(d, e):
            '''
            Расчет общего давления по влагосодержанию и парциальному давлению

            Args:
                d (float или np.ndarray): Влагосодержание [доля]
                e (float или np.ndarray): Парциальное давление водяного пара [Па]

            Returns:
                float или np.ndarray: Атмосферное давление [Па]
            '''
            d = np.asarray(d, dtype=np.float64)
            e = np.asarray(e, dtype=np.float64)
            return e * (M + d) / d

        @staticmethod
        def rho_vT(rho, vT):
            '''
            Расчет общего давления по плотности и виртуальной температуре

            Args:
                rho (float или np.ndarray): Плотность влажного воздуха [кг/м³]
                vT (float или np.ndarray): Виртуальная температура [K]

            Returns:
                float или np.ndarray: Атмосферное давление [Па]
            '''
            rho = np.asarray(rho, dtype=np.float64)
            vT = np.asarray(vT, dtype=np.float64)
            return R_air * vT * rho

    class c:
        '''Методы для работы с удельной теплоемкостью'''

        @staticmethod
        def d(d):
            '''
            Расчет удельной теплоемкости по влагосодержанию

            Args:
                d (float или np.ndarray): Влагосодержание [доля]

            Returns:
                float или np.ndarray: Удельная теплоемкость влажного воздуха [Дж/(кг·K)]
            '''
            d = np.asarray(d, dtype=np.float64)
            return C_p_air * (1 - d) + C_p_h2o * d

    class V:
        '''Методы для работы с удельным объемом'''

        @staticmethod
        def rho(rho):
            '''
            Расчет удельного объема по плотности

            Args:
                rho (float или np.ndarray): Плотность влажного воздуха [кг/м³]

            Returns:
                float или np.ndarray: Удельный объем [м³/кг]
            '''
            rho = np.asarray(rho, dtype=np.float64)
            return 1 / rho


if __name__ == '__main__':
    # Пример расчетов
    print(Humid_air.density(10, 101325, 0.5))  # 1.244425737661975
    print(Humid_air.density([10], np.array(
        [101325]), np.array([0.5])))  # [1.24442581]
    # [1.2467098  1.24613878 1.24556777 1.24499675 1.24442574 1.24385472 1.24328371 1.24271269 1.24214168 1.24157066]
    print(Humid_air.density(10, 101325, np.array(
        [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1])))
