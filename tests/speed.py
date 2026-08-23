import time
import math
import numpy as np


def test_performance():
    """Тест производительности: замеряем инициализацию и вычисления"""

    # 1. Замеряем время импорта модуля
    start_import = time.perf_counter()
    from humid_air import Humid_air
    import_time = time.perf_counter() - start_import
    print(f"⏱ Импорт модуля: {import_time:.4f} сек")

    # 2. Замеряем время первого вызова (компиляция/загрузка кеша)
    start_first = time.perf_counter()
    result1 = Humid_air.density(10, 101325, 0.5)
    first_call_time = time.perf_counter() - start_first
    print(
        f"⏱ Первый вызов (кеширование/компиляция): {first_call_time:.4f} сек")
    print(f"   Результат: {result1}")

    # 3. Замеряем время второго вызова (должно быть быстро)
    start_second = time.perf_counter()
    result2 = Humid_air.density(10, 101325, 0.5)
    second_call_time = time.perf_counter() - start_second
    print(f"⏱ Второй вызов (уже скомпилировано): {second_call_time:.4f} сек")
    print(f"   Результат: {result2}")

    # 4. Тест на массиве из 10 элементов
    data = np.array([0.1, math.nan, 0.3, math.nan, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])

    # Прогрев
    for _ in range(5):
        Humid_air.density(10, 101325, data)

    # Замер 10 повторений
    times = []
    for i in range(10):
        start = time.perf_counter()
        result = Humid_air.density(10, 101325, data)
        times.append(time.perf_counter() - start)

    avg_time = np.mean(times) * 1000  # в миллисекундах
    std_time = np.std(times) * 1000
    print(f"⏱ 10 вычислений на массиве (10 элементов):")
    print(f"   Среднее: {avg_time:.4f} мс")
    print(f"   Std: {std_time:.4f} мс")
    print(f"   Min: {min(times) * 1000:.4f} мс")
    print(f"   Max: {max(times) * 1000:.4f} мс")
    print(
        f"   Результат (первые 3): {result[:3] if isinstance(result, np.ndarray) else result}")

    # 5. Тест на массиве из 1 000 000 элементов
    big_data = np.random.rand(1_000_000) * 0.5 + 0.25  # 0.25-0.75

    # Прогрев
    for _ in range(3):
        Humid_air.density(10, 101325, big_data)

    # Замер 5 повторений
    times_big = []
    for i in range(5):
        start = time.perf_counter()
        result = Humid_air.density(10, 101325, big_data)
        times_big.append(time.perf_counter() - start)

    avg_time_big = np.mean(times_big) * 1000
    std_time_big = np.std(times_big) * 1000
    print(f"⏱ 5 вычислений на массиве (1M элементов):")
    print(f"   Среднее: {avg_time_big:.4f} мс")
    print(f"   Std: {std_time_big:.4f} мс")
    print(f"   Min: {min(times_big) * 1000:.4f} мс")
    print(f"   Max: {max(times_big) * 1000:.4f} мс")

    return {
        'import_time': import_time,
        'first_call': first_call_time,
        'second_call': second_call_time,
        'avg_small': avg_time,
        'avg_big': avg_time_big
    }


if __name__ == "__main__":
    results = test_performance()
