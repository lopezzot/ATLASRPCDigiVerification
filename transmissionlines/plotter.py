import glob
import numpy as np
import matplotlib.pyplot as plt
import re

def extract_time_ns(filename):
    match = re.search(r't([\d\.]+)ns', filename)
    return float(match.group(1)) if match else -1

# Trova tutti i file output_t*.txt
files = sorted(glob.glob("rpc*.txt"), key=extract_time_ns)

plt.figure(figsize=(10, 6))

for f in files:
    time_ns = extract_time_ns(f)
    data = np.loadtxt(f)
    x = data[:, 0]
    V = data[:, 1]
    plt.plot(x, V, label=f'{time_ns:.5f} ns')

plt.xlabel("Strip position (m)")
plt.ylabel("Voltage (V)")
plt.title("RPC signal evolution")
plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.2),
           ncol=3, fontsize='small')  # ncol regolabile
plt.tight_layout()
plt.grid(True)
plt.tight_layout()
plt.show()

