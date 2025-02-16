import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D  # Import for custom legend elements

import mplhep as hep

# Matplotlib configuration for publication-quality plots
hep.style.use("CMS")

# ===============================
#  Read Data from File
# ===============================
file_path = "xs_Table_FM0_FM1.txt"

# Initialize lists
FM0_vals_table = []
FM1_vals_table = []
sigma_vals_table = []
efficiencies_table = []

# Read the file and extract data
with open(file_path, "r") as file:
    lines = file.readlines()[1:]  # Skip header

    for line in lines:
        parts = line.strip().split()  # Split by whitespace
        FM0_vals_table.append(float(parts[0]))
        FM1_vals_table.append(float(parts[1]))
        sigma_vals_table.append(float(parts[2]))   
        efficiencies_table.append(float(parts[3]) / 100.0)   

# Convert lists to numpy arrays
FM0_vals_table = np.array(FM0_vals_table)
FM1_vals_table = np.array(FM1_vals_table)
sigma_vals_table = np.array(sigma_vals_table)
efficiencies_table = np.array(efficiencies_table)

# ===============================
#  Compute σ_limit (Upper Limit)
# ===============================

# Extract SM values (first row in the table)
sigma_SM = sigma_vals_table[0]           
eff_background = efficiencies_table[0]  # Background efficiency

# Compute expected background events
luminosity = 1  # pb^-1
N_bkg = sigma_SM * luminosity * eff_background

# Compute uncertainties
delta_stat = np.sqrt(N_bkg)  # Statistical
delta_sys = 0.10 * N_bkg  # Systematic (10%)
delta_total = np.sqrt(delta_stat**2 + delta_sys**2)  # Total uncertainty

# Compute upper limit on signal events at 95% and 68% CL
N_sig_limit_95 = N_bkg + 1.96 * delta_total  # Gaussian 95% CL
N_sig_limit_68 = N_bkg + 1.00 * delta_total  # Gaussian 68% CL (for one sigma)

# Compute σ_limit using the average signal efficiency
eff_signal_avg = np.mean(efficiencies_table)  # Average signal efficiency
sigma_limit_95 = N_sig_limit_95 / (luminosity * eff_signal_avg)   
sigma_limit_68 = N_sig_limit_68 / (luminosity * eff_signal_avg)   

print(f"🔹 95% CL Upper Limit on Cross-Section: {sigma_limit_95:.6f} pb")
print(f"🔹 68% CL Upper Limit on Cross-Section: {sigma_limit_68:.6f} pb")

# ===============================
#   Define the Exclusion Contour
# ===============================

# Example coefficients from Monte Carlo fits
c0 = 2.001086     # Coefficient for FM0^2
c1 = 0.150223     # Coefficient for FM1^2
c01 = -1.004603   # Mixed interference term

# Constants for transformation
g = 0.628  # Weak coupling constant (example value)
v = 246  # Higgs vacuum expectation value in GeV

# Define the quadratic function for the exclusion limit
def exclusion_contour(FM0, FM1, sigma_limit):
    return c0 * FM0**2 + c1 * FM1**2 + c01 * FM0 * FM1 - (sigma_limit - sigma_SM)

# Generate FM0 and FM1 values
FM0_vals = np.linspace(-3, 3, 200)
FM1_vals = np.linspace(-3, 3, 200)
FM0_grid, FM1_grid = np.meshgrid(FM0_vals, FM1_vals)

# Compute cross-section values for 95% and 68% CL regions
sigma_values_95 = exclusion_contour(FM0_grid, FM1_grid, sigma_limit_95)
sigma_values_68 = exclusion_contour(FM0_grid, FM1_grid, sigma_limit_68)

# Define the energy scale of new physics (TeV)
Lambda = 1  # Ensure this value is correctly defined in your analysis

# Convert FM0 and FM1 to a_0^W and a_C^W using given relations
# Take care of 10^{-8} scaling from the data file
scale_factor = 1e-8  # Apply the correct scaling

a0W_grid = FM0_grid * scale_factor * (Lambda**2 * g**2 * v**2)
aCW_grid = -FM1_grid * scale_factor * (Lambda**2 * g**2 * v**2)


# ===============================
#   Plot the 68% & 95% CL Contours
# ===============================

plt.figure(figsize=(11, 11))
plt.subplots_adjust(left=0.15, right=0.95, bottom=0.12, top=0.95)

# 95% CL Contour (Outer)
contour_95 = plt.contour(a0W_grid, aCW_grid, sigma_values_95, levels=[0], colors='blue', linewidths=2)
plt.contourf(a0W_grid, aCW_grid, sigma_values_95, levels=[-5, 0], colors=['blue'], alpha=0.3)

# 68% CL Contour (Inner)
contour_68 = plt.contour(a0W_grid, aCW_grid, sigma_values_68, levels=[0], colors='magenta', linewidths=2)
plt.contourf(a0W_grid, aCW_grid, sigma_values_68, levels=[-5, 0], colors=['magenta'], alpha=0.5)

# Mark the Standard Model point at (0,0)
plt.scatter(0, 0, color='red', s=30, label="Standard Model")

# Labels and Title
plt.xlabel(r"$a_0^W/\Lambda^2$ [GeV$^{-2}$]")
plt.ylabel(r"$a_C^W/\Lambda^2$ [GeV$^{-2}$]")
plt.title(r"LHeC@1.2 TeV : $e^- p \to e^- w^+ w^- p \to e^- j j \ell \nu_{\ell} p$", pad=15)  # Adjusted title position     (100 fb$^{-1}$)

# Add legend using proper formatting
legend_elements = [
    Line2D([0], [0], color='blue', lw=2, label="95% Confidence Region"),
    Line2D([0], [0], color='magenta', lw=2, label="68% Confidence Region"),
    Line2D([0], [0], marker='o', color='w', markerfacecolor='red', markersize=8, label="Standard Model")
]

plt.legend(handles=legend_elements, loc="upper left", frameon=False, edgecolor='black')

plt.grid()

# Save and show the plot
plt.savefig("a0W_aCW_Exclusion_Contour.jpg", dpi=300)
plt.show()
