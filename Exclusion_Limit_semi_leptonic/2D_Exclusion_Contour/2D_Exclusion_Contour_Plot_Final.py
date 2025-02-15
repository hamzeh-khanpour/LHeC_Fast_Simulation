import numpy as np
import matplotlib.pyplot as plt
import mplhep as hep
from matplotlib.lines import Line2D  # Import for custom legend elements

# Matplotlib configuration for publication-quality plots
hep.style.use("CMS")

# ===============================
#  Read Data from File
# ===============================
file_path = "xs_Table_FM0_FM1.txt"

# Load data (skip header)
data = np.loadtxt(file_path, skiprows=1, dtype=str)

# Extract columns
FM0_vals = data[:, 0].astype(float)  # FM0 values
FM1_vals = data[:, 1].astype(float)  # FM1 values
sigma_vals = data[:, 2].astype(float)  # Cross-section in pb
efficiency_vals = np.array([float(e.strip('%')) / 100 for e in data[:, 3]])  # Convert efficiency from % to decimal

# ===============================
#  Compute σ_limit (Upper Limit)
# ===============================

# Experimental parameters
Lumi = 100  # Integrated Luminosity in fb^-1
Background = 13.35  # Total background events

# Compute S_95 (Upper limit on signal events at 95% CL)
S_95 = 1.96 * np.sqrt(Background)
print(f"S_95 (Upper limit on signal events) = {S_95:.2f}")

# Compute expected signal events
S_vals = sigma_vals * Lumi * efficiency_vals  # Expected signal event count

# Compute upper limit on cross-section using the average signal efficiency
eff_signal_avg = np.mean(efficiency_vals)  
sigma_limit = S_95 / (Lumi * eff_signal_avg)  # Convert to pb

print(f"🔹 95% CL Upper Limit on Cross-Section: {sigma_limit:.6f} pb")

# ===============================
#   Define the Exclusion Contour
# ===============================

# Example quadratic coefficients from Monte Carlo fits
c0 = 2.001086     # Coefficient for FM0^2
c1 = 0.150223     # Coefficient for FM1^2
c01 = -1.004603   # Mixed interference term

# Define the quadratic function for exclusion
def exclusion_contour(FM0, FM1):
    return c0 * FM0**2 + c1 * FM1**2 + c01 * FM0 * FM1 - (sigma_limit - sigma_vals[0])

# Generate FM0 and FM1 values for contour grid
FM0_range = np.linspace(-5, 5, 200)
FM1_range = np.linspace(-5, 5, 200)
FM0_grid, FM1_grid = np.meshgrid(FM0_range, FM1_range)

# Compute exclusion boundary values
sigma_values = exclusion_contour(FM0_grid, FM1_grid)

# ===============================
#   Plot the 95% CL Contour
# ===============================

plt.figure(figsize=(9, 10))
plt.subplots_adjust(left=0.15, right=0.95, bottom=0.12, top=0.95)

# 95% CL Contour
contour_95 = plt.contour(FM0_grid, FM1_grid, sigma_values, levels=[0], colors='blue', linewidths=2)
plt.contourf(FM0_grid, FM1_grid, sigma_values, levels=[-5, 0], colors=['blue'], alpha=0.3)

# Mark the Standard Model point at (0,0)
plt.scatter(0, 0, color='red', s=30, label="Standard Model")

# Labels and Title
plt.xlabel(r"$F_{M0}/\Lambda^4$ [TeV$^{-4}$]")
plt.ylabel(r"$F_{M1}/\Lambda^4$ [TeV$^{-4}$]")
plt.title(r"LHeC @ 1.2 TeV (100 fb$^{-1}$)", pad=15)

# Add legend
legend_elements = [
    Line2D([0], [0], color='blue', lw=2, label="95% Confidence Region"),
    Line2D([0], [0], marker='o', color='w', markerfacecolor='red', markersize=8, label="Standard Model")
]
plt.legend(handles=legend_elements, loc="upper left", frameon=False, edgecolor='black')

plt.grid()

# Save and show the plot
plt.savefig("95CL_Exclusion_Contour_QuadraticFit.jpg", dpi=300)
plt.show()
