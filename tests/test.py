from humid_air_py import Humid_air

print(Humid_air.density(10, 101325, 0.5)) # 1.244425737661975

print(Humid_air.density([10], (101325), [0.5])) # [1.24442581]

print(Humid_air.density(10, 101325, [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]))
# [1.2467098  1.24613878 1.24556777 1.24499675 1.24442574 1.24385472 1.24328371 1.24271269 1.24214168 1.24157066]

t_air = 25
p_air = 101325
G_air = 1
h_air = 0.90
t_water = 15
minT = 5
maxH = 0.95
d = Humid_air.moisture_content(t_air, p_air, h_air)
G_air = G_air + d
v1 = Humid_air.g_water(t_air, p_air, h_air, t_water, minT, maxH)
print(v1)