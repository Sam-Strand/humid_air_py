"""
Тестирование согласованности состояния воздуха.
Для каждой точки сравниваем все возможные пути вычисления параметров.
"""

import json
import pytest
from pathlib import Path
from humid_air import Humid_air

RTOL = 1e-6
ATOL_TEMP = 0.01
ATOL_PRESSURE = 0.1
ATOL_RATIO = 1e-6

def load_snapshots():
    snapshot_path = Path(__file__).parent / "state_snapshots.json"
    if not snapshot_path.exists():
        pytest.fail("Снимки не найдены! Запустите generate_state_snapshots.py")
    
    with open(snapshot_path) as f:
        return json.load(f)

SNAPSHOTS = load_snapshots()

class TestStateConsistency:
    """Проверка, что все пути вычисления дают одинаковые результаты"""
    
    @pytest.mark.parametrize("point_idx", range(len(SNAPSHOTS["points"])))
    def test_E_sat_consistency(self, point_idx):
        """E_sat должно быть одинаковым из разных функций"""
        point = SNAPSHOTS["points"][point_idx]
        state = point["state"]
        inp = state["input"]
        
        # Прямое вычисление
        E_sat_direct = Humid_air.E_sat.t(inp["t"])
        
        # Через e и h
        e = Humid_air.e_evap.E_sat_h(E_sat_direct, inp["h"])
        E_sat_from_e_h = Humid_air.E_sat.e_evap_h(e, inp["h"])
        
        assert E_sat_direct == pytest.approx(E_sat_from_e_h, rel=RTOL, abs=ATOL_PRESSURE)
    
    @pytest.mark.parametrize("point_idx", range(len(SNAPSHOTS["points"])))
    def test_e_evap_consistency(self, point_idx):
        """Парциальное давление должно быть одинаковым из разных функций"""
        point = SNAPSHOTS["points"][point_idx]
        state = point["state"]
        inp = state["input"]
        
        # Базовые значения
        E_sat = Humid_air.E_sat.t(inp["t"])
        d = Humid_air.moisture_content(inp["t"], inp["p"], inp["h"])
        vT = Humid_air.vT.e_p_t(
            Humid_air.e_evap.E_sat_h(E_sat, inp["h"]), 
            inp["p"], inp["t"]
        )
        
        # Три пути вычисления e
        e1 = Humid_air.e_evap.E_sat_h(E_sat, inp["h"])
        e2 = Humid_air.e_evap.d_p(d, inp["p"])
        e3 = Humid_air.e_evap.p_t_vT(inp["p"], inp["t"], vT)
        
        assert e1 == pytest.approx(e2, rel=RTOL, abs=ATOL_PRESSURE)
        assert e1 == pytest.approx(e3, rel=RTOL, abs=ATOL_PRESSURE)
    
    @pytest.mark.parametrize("point_idx", range(len(SNAPSHOTS["points"])))
    def test_d_consistency(self, point_idx):
        """Влагосодержание должно быть одинаковым из разных функций"""
        point = SNAPSHOTS["points"][point_idx]
        state = point["state"]
        inp = state["input"]
        
        d_direct = Humid_air.moisture_content(inp["t"], inp["p"], inp["h"])
        e = Humid_air.e_evap.E_sat_h(Humid_air.E_sat.t(inp["t"]), inp["h"])
        i = Humid_air.i.d_t(d_direct, inp["t"])
        c = Humid_air.c.d(d_direct)
        
        # Разные пути
        d_from_e_p = Humid_air.d.e_evap_p(e, inp["p"])
        d_from_c = Humid_air.d.c(c)
        d_from_i_t = Humid_air.d.i_t(i, inp["t"])
        
        assert d_direct == pytest.approx(d_from_e_p, rel=RTOL, abs=ATOL_RATIO)
        assert d_direct == pytest.approx(d_from_c, rel=RTOL, abs=ATOL_RATIO)
        assert d_direct == pytest.approx(d_from_i_t, rel=RTOL, abs=ATOL_RATIO)
    
    @pytest.mark.parametrize("point_idx", range(len(SNAPSHOTS["points"])))
    def test_t_consistency(self, point_idx):
        """Температура должна быть одинаковой из разных функций"""
        point = SNAPSHOTS["points"][point_idx]
        state = point["state"]
        inp = state["input"]
        
        d = Humid_air.moisture_content(inp["t"], inp["p"], inp["h"])
        e = Humid_air.e_evap.E_sat_h(Humid_air.E_sat.t(inp["t"]), inp["h"])
        E_sat = Humid_air.E_sat.t(inp["t"])
        vT = Humid_air.vT.e_p_t(e, inp["p"], inp["t"])
        i = Humid_air.i.d_t(d, inp["t"])
        L_sat = Humid_air.L_sat.t(inp["t"])
        
        # Разные пути
        t_from_L_sat = Humid_air.t.L_sat(L_sat)
        t_from_e_p_vT = Humid_air.t.e_evap_p_vT(e, inp["p"], vT)
        t_from_E_sat = Humid_air.t.E_sat(E_sat)
        t_from_i_d = Humid_air.t.i_d(i, d)
        
        assert inp["t"] == pytest.approx(t_from_L_sat, rel=RTOL, abs=ATOL_TEMP)
        assert inp["t"] == pytest.approx(t_from_e_p_vT, rel=RTOL, abs=ATOL_TEMP)
        assert inp["t"] == pytest.approx(t_from_E_sat, rel=RTOL, abs=ATOL_TEMP)
        assert inp["t"] == pytest.approx(t_from_i_d, rel=RTOL, abs=ATOL_TEMP)
    
    @pytest.mark.parametrize("point_idx", range(len(SNAPSHOTS["points"])))
    def test_rho_consistency(self, point_idx):
        """Плотность должна быть одинаковой из разных функций"""
        point = SNAPSHOTS["points"][point_idx]
        state = point["state"]
        inp = state["input"]
        
        rho_direct = Humid_air.density(inp["t"], inp["p"], inp["h"])
        e = Humid_air.e_evap.E_sat_h(Humid_air.E_sat.t(inp["t"]), inp["h"])
        vT = Humid_air.vT.e_p_t(e, inp["p"], inp["t"])
        V = 1.0 / rho_direct
        
        rho_from_p_vT = Humid_air.rho.p_vT(inp["p"], vT)
        rho_from_V = Humid_air.rho.V(V)
        
        assert rho_direct == pytest.approx(rho_from_p_vT, rel=RTOL, abs=ATOL_RATIO)
        assert rho_direct == pytest.approx(rho_from_V, rel=RTOL, abs=ATOL_RATIO)
    
    @pytest.mark.parametrize("point_idx", range(len(SNAPSHOTS["points"])))
    def test_vT_consistency(self, point_idx):
        """Виртуальная температура должна быть одинаковой из разных функций"""
        point = SNAPSHOTS["points"][point_idx]
        state = point["state"]
        inp = state["input"]
        
        e = Humid_air.e_evap.E_sat_h(Humid_air.E_sat.t(inp["t"]), inp["h"])
        rho = Humid_air.density(inp["t"], inp["p"], inp["h"])
        
        vT_from_e_p_t = Humid_air.vT.e_p_t(e, inp["p"], inp["t"])
        vT_from_p_rho = Humid_air.vT.p_rho(inp["p"], rho)
        
        assert vT_from_e_p_t == pytest.approx(vT_from_p_rho, rel=RTOL, abs=ATOL_TEMP)
    
    @pytest.mark.parametrize("point_idx", range(len(SNAPSHOTS["points"])))
    def test_p_consistency(self, point_idx):
        """Давление должно быть одинаковым из разных функций"""
        point = SNAPSHOTS["points"][point_idx]
        state = point["state"]
        inp = state["input"]
        
        e = Humid_air.e_evap.E_sat_h(Humid_air.E_sat.t(inp["t"]), inp["h"])
        d = Humid_air.moisture_content(inp["t"], inp["p"], inp["h"])
        vT = Humid_air.vT.e_p_t(e, inp["p"], inp["t"])
        rho = Humid_air.density(inp["t"], inp["p"], inp["h"])
        
        p_from_e_t_vT = Humid_air.p.e_t_vT(e, inp["t"], vT)
        p_from_d_e = Humid_air.p.d_e(d, e)
        p_from_rho_vT = Humid_air.p.rho_vT(rho, vT)
        
        assert inp["p"] == pytest.approx(p_from_e_t_vT, rel=RTOL, abs=ATOL_PRESSURE)
        assert inp["p"] == pytest.approx(p_from_d_e, rel=RTOL, abs=ATOL_PRESSURE)
        assert inp["p"] == pytest.approx(p_from_rho_vT, rel=RTOL, abs=ATOL_PRESSURE)

if __name__ == "__main__":
    # Запуск тестов при прямом вызове
    pytest.main([__file__, "-v", "-s"])
