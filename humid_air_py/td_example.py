import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter

import numpy as np
from humid_air_py import Humid_air


def create_td_diagram(p=101325, t_min=-10, t_max=50, d_max=0.04, figsize=(12, 8)):
    '''
    Построение t-d диаграммы влажного воздуха
    '''
    fig, ax = plt.subplots(figsize=figsize)
    
    # Расчет максимального влагосодержания
    E = Humid_air.E.t(t_max)
    max_d = Humid_air.d.e_p(E, p)
    d_max = min(d_max, max_d)
    
    # Генерация сетки
    t_range = np.linspace(t_min, t_max, 100)
    d_range = np.linspace(0, d_max, 100)
    
    # 1. Линии постоянной относительной влажности
    humidity_levels = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
    colors_h = plt.cm.Blues(np.linspace(0.3, 1, len(humidity_levels)))
    
    for h, color in zip(humidity_levels, colors_h):
        d_values = []
        valid_t = []
        
        for t in t_range:
            try:
                E = Humid_air.E.t(t)
                e = Humid_air.e.E_h(E, h)
                d = Humid_air.d.e_p(e, p)
                if d <= d_max:
                    d_values.append(d)
                    valid_t.append(t)
            except:
                continue
        
        if d_values:
            ax.plot(np.array(d_values) * 1000, valid_t, color=color, linewidth=1.5, label=f'φ={int(h*100)}%')
    
    # 2. Линии постоянной температуры
    t_const_levels = range(int(t_min), int(t_max)+1, 5)
    for t_val in t_const_levels:
        if t_val >= t_min and t_val <= t_max:
            ax.plot(d_range * 1000, [t_val] * len(d_range), 'r--', alpha=0.5, linewidth=0.8)
    
    # 3. Линии постоянного влагосодержания
    d_const_levels = np.arange(0.005, d_max, 0.005)
    for d_val in d_const_levels:
        ax.plot([d_val * 1000] * len(t_range), t_range, 'g--', alpha=0.5, linewidth=0.8)
    
    # 4. Линии постоянной энтальпии
    i_min = Humid_air.i.d_t(0, t_min)
    i_max = Humid_air.i.d_t(d_max, t_max)
    i_levels = np.arange(i_min // 10000 * 10000, i_max, 20000)
    
    for i_val in i_levels:
        d_values = []
        valid_t = []
        
        # Итерируем от высокой температуры к низкой
        for t in np.linspace(t_max, t_min, 100):
            try:
                d = Humid_air.d.i_t(i_val, t)
                if 0 <= d <= d_max:
                    d_values.append(d)
                    valid_t.append(t)
            except:
                continue
        
        if d_values:
            ax.plot(np.array(d_values) * 1000, valid_t, 'purple', alpha=0.7, linewidth=1)
            ax.text(d_values[-1] * 1000, valid_t[-1], f'{i_val/1000:.0f} кДж', 
                       fontsize=8, color='purple', ha='left', va='center')
    
    ax.set_xlabel('Влагосодержание d, г/кг', fontsize=12)
    ax.set_ylabel('Температура t, °C', fontsize=12)
    ax.set_title(f't-d диаграмма влажного воздуха (p = {p/100:.0f} гПа)', 
                fontsize=14, fontweight='bold')
    
    ax.grid(True, alpha=0.3)
    ax.legend(loc='upper left', bbox_to_anchor=(1.05, 1))
    ax.set_xlim(0, d_max * 1000)
    ax.set_ylim(t_min, t_max)

    ax.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))
    
    plt.tight_layout()
    return fig, ax

if __name__ == '__main__':
    fig, ax = create_td_diagram()
    plt.show()