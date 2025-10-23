import matplotlib.pyplot as plt
import glob
import os

# Trova tutti i file che iniziano con "secondcoordinate" e terminano con ".txt"
files = sorted(glob.glob("secondcoordinate*.txt"))

for filename in files:
    x_meas_values = []
    x_real = None
    jitter = None
    tdc_size = None

    with open(filename, 'r') as f:
        for line in f:
            if line.strip() == "":
                continue
            if line.startswith("X_real"):
                continue  # salta l'intestazione
            parts = line.strip().split()
            if len(parts) != 4:
                continue  # salta righe mal formattate

            x_r, x_m, jit, tdc = map(float, parts)

            if x_real is None:
                x_real = x_r
                jitter = jit
                tdc_size = tdc

            x_meas_values.append(x_m)

    # Calcolo dei limiti estesi
    min_x = min(x_meas_values)
    max_x = max(x_meas_values)
    margin = 0.1 * (max_x - min_x)
    xlim_low = min_x - margin
    xlim_high = max_x + margin

    # Crea il plot per questo file
    plt.figure(figsize=(10, 6))
    plt.hist(x_meas_values, bins=30, color='skyblue', edgecolor='black')

    label = f"Pos = {x_real:.3f} m\njitter = {jitter:.3f} ns\nTDC size = {tdc_size:.3f} ns"
    plt.legend([label], loc="upper right")

    plt.xlabel("Second coordinate measurement (m)")
    plt.ylabel("Frequency")
    plt.title(f"")
    plt.xlim(xlim_low, xlim_high)
    plt.grid(True)
    plt.tight_layout()

    plt.show()

