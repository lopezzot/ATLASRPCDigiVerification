import numpy as np
import matplotlib.pyplot as plt

# Carica direttamente il file saltando la prima riga
data = np.loadtxt("threshold.txt", skiprows=1)

x = data[:, 4]
V = data[:, 7]
Y = data[:, 5]

plt.figure(figsize=(10, 6))
plt.plot(x, V)

plt.xlabel("Threshold (V)")
plt.ylabel("Time over threshold (ns) [left boundary]")
plt.title("Time over threshold as a function of the threshold")
plt.grid(True)
plt.tight_layout()
plt.show()

plt.figure(figsize=(10, 6))
plt.plot(x, Y)

plt.xlabel("Threshold (V)")
plt.ylabel("Time of arrival (ns) [left boundary]")
plt.title("Time of arrival as a function of the threshold")
plt.grid(True)
plt.tight_layout()
plt.show()

