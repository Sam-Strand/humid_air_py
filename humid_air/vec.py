from numbers import Real
from numpy.typing import NDArray
import numpy as np
from typing import TypeVar
from numba import njit

Number = TypeVar('Number', Real, NDArray[np.floating])

def vec(cache=True):
    return njit(cache=cache)
