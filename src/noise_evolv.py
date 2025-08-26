# used for 7/23/25 chamber commissioning slides
# but in case somehow it is needed in the future
# Cinyu Zhu  7/23/25

import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 13})

# Dates and x-axis positions
dates = ['July 3 \n (start)', 'July 10 \n (acopian power supply)', 'July 17 \n (passive balancer)', 'July 22 \n (gen1AVC, current)']
x = range(len(dates))

# PD06 data: fully available
pd06 = {
    'CH0 (SR)': [10.48, 9.62, 6, 5.39],
    'CH1': [7.36, 6.60, 5.3, 4.64],
    'CH2': [9.17, 7.42, 5.8, 4.83],
    'CH3': [7.78, 7.01, 5.2, 5.06]
}

# PD07 data: July 9 missing
pd07 = {
    'CH0': [38.97, None, None,4.63],
    'CH1': [24.69, None, None,4.83],
    'CH2': [27.32, None, None,4.50],
    'CH3': [35.16, None, None,5.09]
}

def plot_noise(data, title, filename):
    plt.figure(figsize=(8, 5))
    for i, (ch, y_vals) in enumerate(data.items()):
        # Build x and y without missing data
        x_clean = [j for j, val in enumerate(y_vals) if val is not None]
        y_clean = [val for val in y_vals if val is not None]
        plt.plot(x_clean, y_clean, marker='o', label=ch, color=f'C{i}')
    plt.xticks(range(len(dates)), dates)
    plt.ylabel("Noise (ADU)")
    plt.title(title)
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(filename, dpi=300)
    plt.show()

# Plot PD06 and PD07
plot_noise(pd06, "PD06: Noise vs Date", "PD06.png")
plot_noise(pd07, "PD07: Noise vs Date", "PD07.png")
