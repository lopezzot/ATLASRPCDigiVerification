import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# Carica i dati
data = np.loadtxt("jpeak.txt", skiprows=1)

x = data[:, 9]
V = data[:, 7]  # TOT
Y = data[:, 5]  # TOA

# Calcola spike charge fraction (normalizzato al massimo Jpeak = 0.007 A/m)
absx = np.abs(x) / 0.007

# Funzione di fit: polinomio di secondo grado
def poly2(x, a, b, c):
    return a * x**2 + b * x + c

# Fit per TOT
popt_V, _ = curve_fit(poly2, absx, V)
fit_V = poly2(absx, *popt_V)
# Generate formula as string
formula_V = f"f(x) = {popt_V[0]:+.3f}·x² {popt_V[1]:+.3f}·x {popt_V[2]:+.3f}"

# Fit per TOA
popt_Y, _ = curve_fit(poly2, absx, Y)
fit_Y = poly2(absx, *popt_Y)

# --- Plot TOT ---
plt.figure(figsize=(10, 6))
plt.plot(absx, V, 'o', markersize=4, label="Data TOT")
plt.plot(absx, fit_V, '-', color='red', label="Fit TOT")
plt.xlabel("Spike charge fraction")
plt.ylabel("Time over threshold (ns) [left boundary]")
plt.title("TOT vs. spike charge fraction")
plt.grid(True)
plt.legend()
plt.legend()
plt.text(0.25, 0.95, formula_V, transform=plt.gca().transAxes,
         fontsize=10, verticalalignment='top', bbox=dict(facecolor='white', alpha=0.7))
plt.tight_layout()
plt.show()

# --- Plot TOA ---
plt.figure(figsize=(10, 6))
plt.plot(absx, Y, 'o', markersize=4, label="Data TOA")
plt.plot(absx, fit_Y, '-', color='green', label="Fit TOA")
plt.xlabel("Spike charge fraction")
plt.ylabel("Time of arrival (ns) [left boundary]")
plt.title("TOA vs. spike charge fraction")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()

