import numpy as np
import matplotlib.pyplot as plt

# Carica direttamente il file saltando la prima riga
data = np.loadtxt("scope_output.txt", skiprows=1)

x = data[:, 0]
Vleft = data[:, 1]
Vright = data[:, 2]

plt.figure(figsize=(10, 6))
plt.plot(x, Vleft)

plt.xlabel("Time (ns)")
plt.ylabel("Voltage left-side (V)")
plt.title("Bare signal at scope (left side)")
plt.grid(True)
plt.tight_layout()
plt.show()

plt.figure(figsize=(10, 6))
plt.plot(x, Vright)

plt.xlabel("Time (ns)")
plt.ylabel("Voltage right-side (V)")
plt.title("Bare signal at scope (right side)")
plt.grid(True)
plt.tight_layout()
plt.show()
