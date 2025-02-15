import numpy as np
import matplotlib.pyplot as plt

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

# Compute upper limit on signal events at 95% CL
N_sig_limit = N_bkg + 1.96 * delta_total  # Gaussian 95% CL

# Compute σ_limit using the average signal efficiency
eff_signal_avg = np.mean(efficiencies_table)  # Average signal efficiency
sigma_limit = N_sig_limit / (luminosity * eff_signal_avg)  # Convert to pb

print(f"🔹 95% CL Upper Limit on Cross-Section: {sigma_limit:.6f} pb")


# Convert σ_limit to fb for consistency with contour calculations
#sigma_limit *= 1000.0  # Convert pb → fb
#sigma_SM *= 1000.0  # Convert pb → fb

# ===============================
#   Define the Exclusion Contour
# ===============================

# Example coefficients from Monte Carlo fits
c0 = 2.001086     # Coefficient for FM0^2
c1 = 0.150223     # Coefficient for FM1^2
c01 = -1.004603   # Mixed interference term

# Define the quadratic function for the exclusion limit
def exclusion_contour(FM0, FM1):
    return c0 * FM0**2 + c1 * FM1**2 + c01 * FM0 * FM1 - (sigma_limit - sigma_SM)

# Generate FM0 and FM1 values
FM0_vals = np.linspace(-3, 3, 100)
FM1_vals = np.linspace(-3, 3, 100)
FM0_grid, FM1_grid = np.meshgrid(FM0_vals, FM1_vals)

# Compute cross-section values
sigma_values = exclusion_contour(FM0_grid, FM1_grid)



# ===============================
#   Plot the 95% CL Contour
# ===============================

plt.figure(figsize=(9, 10))
plt.subplots_adjust(left=0.15, right=0.95, bottom=0.12, top=0.95)

contour = plt.contour(FM0_grid, FM1_grid, sigma_values, levels=[0], colors='green', linewidths=2)
#plt.contourf(FM0_grid, FM1_grid, sigma_values, levels=[-5, 0, 5], colors=['yellow', 'green', 'white'], alpha=0.5)

# Fill inside the contour with yellow
plt.contourf(FM0_grid, FM1_grid, sigma_values, levels=[-5, 0], colors=['yellow'], alpha=0.5)


# Mark the Standard Model point at (0,0)
plt.scatter(0, 0, color='red', s=30, label="Standard Model")

# Labels and Title
plt.xlabel(r"$F_{M0}/\Lambda^4$ [TeV$^{-4}$]")
plt.ylabel(r"$F_{M1}/\Lambda^4$ [TeV$^{-4}$]")
plt.title(r"LHeC @ 1.2 TeV")

# Add legend
plt.legend(["95% Confidence Region"], loc="upper right")


plt.grid()

# Save and show the plot
plt.savefig("95CL_Exclusion_Contour.jpg", dpi=300)
plt.show()
