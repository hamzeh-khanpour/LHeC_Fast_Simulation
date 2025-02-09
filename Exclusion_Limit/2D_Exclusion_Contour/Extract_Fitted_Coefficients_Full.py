import numpy as np
from scipy.optimize import curve_fit

# Load data from text file
data = np.loadtxt("xs_Table_FM0_FM1.txt", skiprows=1)  # Skip header row

# Extract first row separately (Standard Model cross-section)
sigma_SM_value = data[0, 2]  # Third column, first row (σ_SM)

# Remove the first row from data before fitting
data = data[1:]  # Keep only the EFT-dependent rows

# Extract columns for fitting
FM0_vals = data[:, 0]  # First column: FM0 values
FM1_vals = data[:, 1]  # Second column: FM1 values
sigma_vals = data[:, 2]  # Third column: Cross-section values

# Define quadratic function for cross-section (using fixed σ_SM_value)
def cross_section_model(X, c0, c1, c01):
    FM0, FM1 = X  # Unpack input
    return sigma_SM_value + c0 * FM0**2 + c1 * FM1**2 + c01 * FM0 * FM1

# Initial parameter guess (excluding σ_SM)
initial_guess = [1, 1, 1]

# Perform curve fitting (excluding σ_SM)
params, covariance = curve_fit(cross_section_model, (FM0_vals, FM1_vals), sigma_vals, p0=initial_guess)

# Extract fitted parameters
c0_fit, c1_fit, c01_fit = params
param_errors = np.sqrt(np.diag(covariance))  # Uncertainty estimation

# Compute residuals and chi-squared per degree of freedom
sigma_pred = cross_section_model((FM0_vals, FM1_vals), *params)
residuals = sigma_vals - sigma_pred
chi_squared = np.sum((residuals**2) / sigma_pred)  # Approximate uncertainties
dof = len(sigma_vals) - len(params)  # Degrees of freedom
chi_squared_per_dof = chi_squared / dof

# Print results with formatting
print("\n=== Fit Results ===")
print(f"c0 = {c0_fit:.6f} ± {param_errors[0]:.6f}")
print(f"c1 = {c1_fit:.6f} ± {param_errors[1]:.6f}")
print(f"c01 = {c01_fit:.6f} ± {param_errors[2]:.6f}")
print(f"sigma_SM (Fixed) = {sigma_SM_value:.6f} (from input)")
print(f"Chi-squared per degree of freedom: {chi_squared_per_dof:.3f}")
