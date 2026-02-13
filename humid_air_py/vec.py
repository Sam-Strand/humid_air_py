from numba import vectorize


def vec(n, nopython=True, cache=True):
    '''n аргументов float64, возвращает float64'''
    return vectorize([f"float64({', '.join(['float64'] * n)})"], nopython=nopython, cache=cache)
