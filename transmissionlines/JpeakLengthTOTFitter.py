import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
from mpl_toolkits.mplot3d import Axes3D  # Needed to enable 3D projection
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import PolynomialFeatures

# Carica i dati
data = np.loadtxt("jpeak_length.txt", skiprows=1)

J = np.abs(data[:, 9])   # Spike current density (A/m)
J = [x/0.007 for x in J] # move to fraction of spike charge on strip
TOT = data[:, 7]         # TOT left (ns)
D = data[:, 0] / 2       # Half distance to edge (m)

# Griglia regolare
xi = np.linspace(min(D), max(D), 200)
yi = np.linspace(min(J), max(J), 200)
X, Y = np.meshgrid(xi, yi)

# Interpolazione TOT su griglia
Z = griddata((D, J), TOT, (X, Y), method='linear')

# Crea figura 3D
fig = plt.figure(figsize=(10, 7))
ax = fig.add_subplot(111, projection='3d')

# Plot della superficie
surf = ax.plot_surface(X, Y, Z, cmap='viridis', edgecolor='none')

# Etichette e colorbar
ax.set_xlabel("Distance to edge (m)")
ax.set_ylabel("Fraction of spike charge on strip")
ax.set_zlabel("Time over threshold (ns)")
ax.set_title("")

fig.colorbar(surf, ax=ax, shrink=0.6, aspect=10, label="")

plt.tight_layout()
#plt.show()

fitD = (data[:, 0] / 2).reshape(-1, 1)   # Distance
fitJ = np.array(J).reshape(-1, 1) # Spike current density
fitTOT = data[:, 7]                      # Time over threshold
# Combine input variables
X_input = np.hstack((fitD, fitJ))
# Generate polynomial features (degree 2 example)
poly = PolynomialFeatures(degree=2, include_bias=False)
X_poly = poly.fit_transform(X_input)
# Fit linear regression on polynomial features
model = LinearRegression()
model.fit(X_poly, fitTOT)
# Calcola TOT predetti sulla stessa griglia
XY_grid = np.stack([X.ravel(), Y.ravel()], axis=1)
XY_poly = poly.transform(XY_grid)
Z_fit = model.predict(XY_poly).reshape(X.shape)
# Output coefficients
coefs = model.coef_
intercept = model.intercept_
feature_names = poly.get_feature_names_out(["D", "J"])
# Display formula
for name, coef in zip(feature_names, coefs):
    print(f"{coef:+.4f} * {name}")
print(f"Intercept: {intercept:+.4f}")

#ax.plot_surface(X, Y, Z_fit, cmap='plasma', alpha=0.6, edgecolor='none') #overlay fit to data
plt.show()

# plot the difference w.r.t. the fit
plt.figure(figsize=(10, 7))
diff = Z - Z_fit
diff = np.nan_to_num(diff, nan=0.0)

plt.title("Difference between interpolated data and fitted model")
plt.xlabel("Distance to edge (m)")
plt.ylabel("Fraction of spike charge on strip")

# Usa pcolormesh per la differenza sulla stessa griglia 2D
pcm = plt.pcolormesh(X, Y, diff, cmap='coolwarm', shading='auto')
plt.colorbar(pcm, label="Difference (ns)")

plt.grid(True)
plt.tight_layout()
plt.show()
