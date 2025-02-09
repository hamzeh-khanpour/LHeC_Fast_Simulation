import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# Load data from text file (skip header row)
data = np.loadtxt("xs_Table_FM0_FM1.txt", skiprows=1, dtype=str)

# Extract columns
FM0_vals = data[:, 0].astype(float)  # FM0 values
FM1_vals = data[:, 1].astype(float)  # FM1 values
sigma_vals = data[:, 2].astype(float)  # Cross-section in pb
efficiency_vals = np.array([float(e.strip('%')) / 100 for e in data[:, 3]])  # Convert efficiency from % to decimal

# Experimental parameters
Lumi = 100  # Integrated Luminosity in fb^-1
Background = 13.35  # Total number of background events

# Compute S_95 (95% Confidence Level limit on signal events)
S_95 = 1.96 * np.sqrt(Background)
print(f"S_95 (Upper limit on signal events) = {S_95:.2f}")

# Compute expected signal events for each (FM0, FM1)
S_vals = sigma_vals * Lumi * efficiency_vals  # Signal event count

# Determine exclusion condition: Exclude if S > S_95
excluded = S_vals > S_95

# Reshape data for contour plotting
FM0_grid = np.unique(FM0_vals)
FM1_grid = np.unique(FM1_vals)
FM0_mesh, FM1_mesh = np.meshgrid(FM0_grid, FM1_grid)

# Create exclusion matrix
exclusion_matrix = np.zeros_like(FM0_mesh, dtype=int)

# Populate exclusion matrix
for i in range(len(FM0_vals)):
    f0_idx = np.where(FM0_grid == FM0_vals[i])[0][0]
    f1_idx = np.where(FM1_grid == FM1_vals[i])[0][0]
    exclusion_matrix[f1_idx, f0_idx] = excluded[i]

# Generate contour plot
plt.figure(figsize=(8,6))

plt.contourf(FM0_mesh, FM1_mesh, exclusion_matrix, levels=[-0.1, 0.1, 1.1], colors=['yellow', 'green'], alpha=0.6)
contour = plt.contour(FM0_mesh, FM1_mesh, exclusion_matrix, levels=[0.5], colors='black', linewidths=2)

# Labels and title
plt.xlabel(r"$F_{M0}/\Lambda^4$ [$\text{TeV}^{-4}$]")
plt.ylabel(r"$F_{M1}/\Lambda^4$ [$\text{TeV}^{-4}$]")
plt.title(r"95% CL Exclusion Contour for $F_{M0}$ and $F_{M1}$")
plt.grid()
plt.show()
