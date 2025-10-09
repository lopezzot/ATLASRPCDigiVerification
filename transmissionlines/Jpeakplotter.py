import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# Carica i dati
data = np.loadtxt("jpeak.txt", skiprows=1)

x = data[:, 9]  # Charge on strip (fC)
V = data[:, 7]  # TOT
Y = data[:, 5]  # TOA

# Funzione di fit: polinomio di secondo grado
def poly2(x, a, b, c):
    return a * x**2 + b * x + c

# Fit per TOT
popt_V, _ = curve_fit(poly2, x, V)
fit_V = poly2(x, *popt_V)
# Generate formula as string
formula_V = f"f(x) = {popt_V[0]:+.5f}·x² {popt_V[1]:+.5f}·x {popt_V[2]:+.5f}"

# Fit per TOA
popt_Y, _ = curve_fit(poly2, x, Y)
fit_Y = poly2(x, *popt_Y)

# --- Plot TOT ---
plt.figure(figsize=(10, 6))
plt.plot(x, V, 'o', markersize=4, label="Data TOT")
plt.plot(x, fit_V, '-', color='red', label="Fit TOT")
plt.xlabel("Charge on strip (fC)")
plt.ylabel("Time over threshold (ns) [left boundary]")
plt.title("TOT vs. charge on strip")
plt.grid(True)
plt.legend()
plt.legend()
plt.text(0.25, 0.95, formula_V, transform=plt.gca().transAxes,
         fontsize=10, verticalalignment='top', bbox=dict(facecolor='white', alpha=0.7))
plt.tight_layout()
plt.show()

# --- Plot TOA ---
plt.figure(figsize=(10, 6))
plt.plot(x, Y, 'o', markersize=4, label="Data TOA")
plt.plot(x, fit_Y, '-', color='green', label="Fit TOA")
plt.xlabel("Charge on strip (fC)")
plt.ylabel("Time of arrival (ns) [left boundary]")
plt.title("TOA vs. charge on strip")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()

