import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

# Define the quadratic function for the exclusion limit
def exclusion_contour(FM0, FM1, c0, c1, c01, sigma_limit, sigma_SM):
    return c0 * FM0**2 + c1 * FM1**2 + c01 * FM0 * FM1 - (sigma_limit - sigma_SM)

# Example coefficients (these should be obtained from actual MC fits)
c0 = 2.001086  # Coefficient for FM0^2
c1 = 0.150223  # Coefficient for FM1^2
c01 = -1.004603  # Mixed interference term

sigma_SM = 0.0134936 * 1000.0  # SM cross-section in fb
sigma_limit = 0.0154936  * 1000.0  # Upper limit on the cross-section

# Generate FM0 and FM1 values
FM0_vals = np.linspace(-20, 20, 100)
FM1_vals = np.linspace(-20, 20, 100)
FM0_grid, FM1_grid = np.meshgrid(FM0_vals, FM1_vals)

# Compute cross-section values
sigma_values = exclusion_contour(FM0_grid, FM1_grid, c0, c1, c01, sigma_limit, sigma_SM)

# Create the contour plot
plt.figure(figsize=(8,6))
contour = plt.contour(FM0_grid, FM1_grid, sigma_values, levels=[0], colors='black', linewidths=2)
plt.contourf(FM0_grid, FM1_grid, sigma_values, levels=[-20, 0, 20], colors=['yellow', 'green', 'white'], alpha=0.5)

# Labels and Title
plt.xlabel(r"$F_{M0}/\Lambda^4$ [$\text{TeV}^{-4}$]")
plt.ylabel(r"$F_{M1}/\Lambda^4$ [$\text{TeV}^{-4}$]")
plt.title(r"95% CL Exclusion Contour on $F_{M0}$ and $F_{M1}$")
plt.grid()
plt.show()
