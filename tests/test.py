import math
from humid_air import Humid_air
import numpy as np

print(Humid_air.density(10, 101325, 0.5)) # 1.244425737661975

print(Humid_air.density(np.array([10]), np.array((101325)), np.array([0.5]))) # [1.24442581]

print(Humid_air.density(10, 101325, np.array([0.1, math.nan, 0.3, math.nan, 0.5, 0.6, 0.7, 0.8, 0.9, 1])))
# [1.2467098  1.24613878 1.24556777 1.24499675 1.24442574 1.24385472 1.24328371 1.24271269 1.24214168 1.24157066]
