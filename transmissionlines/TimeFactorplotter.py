import numpy as np
import matplotlib.pyplot as plt

# Carica direttamente il file saltando la prima riga
data = np.loadtxt("timefactor.txt", skiprows=1)

x = data[:, 0]
V = data[:, 1]

plt.figure(figsize=(10, 6))
plt.plot(x, V)

plt.xlabel("Time (ns)")
plt.ylabel("Current density (a.u.)")
plt.title("Time evolution of RPC spike current density")
plt.grid(True)
plt.tight_layout()
plt.show()

