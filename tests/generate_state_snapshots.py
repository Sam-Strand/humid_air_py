"""
Генерация эталонных состояний влажного воздуха.
Создает 100 точек с разными t, p, h и вычисляет ВСЕ параметры
через разные функции для проверки согласованности.
"""

import json
import numpy as np
from humid_air import Humid_air

# Константы точности
RTOL = 1e-6
ATOL_TEMP = 0.01
ATOL_PRESSURE = 0.1
ATOL_RATIO = 1e-6

class NumpyEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        if isinstance(obj, np.floating):
            return float(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return super().default(obj)

def generate_test_points(n_points=100):
    """
    Генерирует n_points равномерно распределенных комбинаций
    t, p, h в допустимых диапазонах
    """
    np.random.seed(42)  # для воспроизводимости
    
    # Диапазоны
    t_range = (-30, 40)      # температура °C
    p_range = (90000, 110000) # давление Па
    h_range = (0.1, 0.9)      # влажность (избегаем 0 и 1)
    
    points = []
    for _ in range(n_points):
        t = np.random.uniform(*t_range)
        p = np.random.uniform(*p_range)
        h = np.random.uniform(*h_range)
        points.append({
            "t": float(t),
            "p": float(p),
            "h": float(h),
            "name": f"point_{_:03d}_t{t:.1f}_p{p:.0f}_h{h:.2f}"
        })
    
    return points

def compute_state(t, p, h):
    """
    Вычисляет ПОЛНОЕ состояние воздуха через разные функции
    Возвращает словарь со всеми параметрами, вычисленными разными способами
    """
    state = {
        "input": {"t": t, "p": p, "h": h},
        
        # Прямые вычисления (базовые)
        "direct": {
            "d": Humid_air.moisture_content(t, p, h),
            "rho": Humid_air.density(t, p, h),
            "t_dew": Humid_air.dew_point_temperature(t, h),
            "d_max": Humid_air.maximum_moisture_content(p, t),
        },
        
        # E_sat (давление насыщения) через разные функции
        "E_sat": {
            "from_t": float(Humid_air.E_sat.t(t)),  # E_sat(t)
            "from_e_h": None,  # будет заполнено позже
        },
        
        # e_evap (парциальное давление) через разные функции
        "e_evap": {
            "from_E_sat_h": None,  # e = e_evap(E_sat, h)
            "from_d_p": None,       # e = e_evap(d, p)
            "from_p_t_vT": None,    # e = e_evap(p, t, vT)
        },
        
        # d (влагосодержание) через разные функции
        "d": {
            "from_e_p": None,  # d = d(e, p)
            "from_L_sat_l_evap": None,  # d = d(L_sat, l_evap)
            "from_c": None,    # d = d(c)
            "from_i_t": None,  # d = d(i, t)
        },
        
        # t (температура) через разные функции
        "t": {
            "from_L_sat": None,        # t = t(L_sat)
            "from_e_p_vT": None,   # t = t(e, p, vT)
            "from_E_sat": None,     # t = t(E_sat)
            "from_i_d": None,       # t = t(i, d)
        },
        
        # rho (плотность) через разные функции
        "rho": {
            "from_p_vT": None,  # rho = rho(p, vT)
            "from_V": None,     # rho = rho(V)
        },
        
        # h (влажность) через разные функции
        "h": {
            "from_E_sat_e": None,  # h = h(E_sat, e)
        },
        
        # i (энтальпия) через разные функции
        "i": {
            "from_d_t": None,  # i = i(d, t)
        },
        
        # vT (виртуальная температура) через разные функции
        "vT": {
            "from_e_p_t": None,  # vT = vT(e, p, t)
            "from_p_rho": None,   # vT = vT(p, rho)
        },
        
        # p (давление) через разные функции
        "p": {
            "from_e_t_vT": None,  # p = p(e, t, vT)
            "from_d_e": None,      # p = p(d, e)
            "from_rho_vT": None,   # p = p(rho, vT)
        },
        
        # c (теплоемкость) через разные функции
        "c": {
            "from_d": None,  # c = c(d)
        },
        
        # V (объем) через разные функции
        "V": {
            "from_rho": None,  # V = V(rho)
        },
    }
    
    # Получаем базовые значения
    d = state["direct"]["d"]
    E_sat = state["E_sat"]["from_t"]
    e = Humid_air.e_evap.E_sat_h(E_sat, h)
    
    # Виртуальная температура
    vT = Humid_air.vT.e_p_t(e, p, t)
    
    # Энтальпия
    i = Humid_air.i.d_t(d, t)
    
    # Теплоемкость
    c = Humid_air.c.d(d)
    
    # Объем для 1 кг воздуха
    V = 1.0 / state["direct"]["rho"]
    
    # Теперь заполняем все перекрестные вычисления
    
    # E_sat
    state["E_sat"]["from_e_h"] = float(Humid_air.E_sat.e_evap_h(e, h))
    
    # e_evap
    state["e_evap"]["from_E_sat_h"] = float(Humid_air.e_evap.E_sat_h(E_sat, h))
    state["e_evap"]["from_d_p"] = float(Humid_air.e_evap.d_p(d, p))
    state["e_evap"]["from_p_t_vT"] = float(Humid_air.e_evap.p_t_vT(p, t, vT))
    
    # d
    state["d"]["from_e_p"] = float(Humid_air.d.e_evap_p(e, p))
    # L_sat и l_evap нужны для d через L_l
    L_sat = Humid_air.L_sat.t(t)
    l_evap = Humid_air.l_evap.d_L_sat(d, L_sat) if d > 0 else 0
    state["d"]["from_L_sat_l_evap"] = float(Humid_air.d.L_sat_l_evap(L_sat, l_evap))
    state["d"]["from_c"] = float(Humid_air.d.c(c))
    state["d"]["from_i_t"] = float(Humid_air.d.i_t(i, t))
    
    # t
    state["t"]["from_L_sat"] = float(Humid_air.t.L_sat(L_sat))
    state["t"]["from_e_p_vT"] = float(Humid_air.t.e_evap_p_vT(e, p, vT))
    state["t"]["from_E_sat"] = float(Humid_air.t.E_sat(E_sat))
    state["t"]["from_i_d"] = float(Humid_air.t.i_d(i, d))
    
    # rho
    state["rho"]["from_p_vT"] = float(Humid_air.rho.p_vT(p, vT))
    state["rho"]["from_V"] = float(Humid_air.rho.V(V))
    
    # h
    state["h"]["from_E_sat_e"] = float(Humid_air.h.E_sat_e_evap(E_sat, e))
    
    # i
    state["i"]["from_d_t"] = float(Humid_air.i.d_t(d, t))
    
    # vT
    state["vT"]["from_e_p_t"] = float(Humid_air.vT.e_p_t(e, p, t))
    state["vT"]["from_p_rho"] = float(Humid_air.vT.p_rho(p, state["direct"]["rho"]))
    
    # p
    state["p"]["from_e_t_vT"] = float(Humid_air.p.e_t_vT(e, t, vT))
    state["p"]["from_d_e"] = float(Humid_air.p.d_e(d, e))
    state["p"]["from_rho_vT"] = float(Humid_air.p.rho_vT(state["direct"]["rho"], vT))
    
    # c
    state["c"]["from_d"] = float(Humid_air.c.d(d))
    
    # V
    state["V"]["from_rho"] = float(Humid_air.V.rho(state["direct"]["rho"]))
    
    return state

def main():
    print("Генерация 100 тестовых точек...")
    points = generate_test_points(20)
    
    snapshots = {
        "metadata": {
            "description": "Полные состояния влажного воздуха для проверки согласованности",
            "n_points": len(points),
            "rtol": RTOL,
            "atol_temp": ATOL_TEMP,
            "atol_pressure": ATOL_PRESSURE,
            "atol_ratio": ATOL_RATIO,
        },
        "points": []
    }
    
    for i, point in enumerate(points):
        print(f"Точка {i+1}/100: t={point['t']:.1f}°C, p={point['p']:.0f}Па, h={point['h']:.2f}")
        state = compute_state(point["t"], point["p"], point["h"])
        
        # Добавляем имя точки
        point_with_state = {
            "name": point["name"],
            "state": state
        }
        snapshots["points"].append(point_with_state)
    
    # Сохраняем
    output_file = "tests/state_snapshots.json"
    with open(output_file, "w") as f:
        json.dump(snapshots, f, indent=2, cls=NumpyEncoder)
    
    print(f"\n✅ {len(points)} состояний сохранено в {output_file}")

if __name__ == "__main__":
    main()