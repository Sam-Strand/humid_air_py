from humid_air import Humid_air as ha

# Расчет точки росы
t_dry = 25.0  # °C, температура сухого термометра
rh = 0.6      # доля, относительная влажность (0.0-1.0)

t_dew = ha.dew_point_temperature(t_dry, rh)
print(f"Точка росы: {t_dew:.2f} °C")
# Точка росы: 16.69 °C