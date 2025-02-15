import numpy as np
from scipy.optimize import curve_fit

# Monte Carlo results (Example data - Replace with real MC values)
FM0_vals = np.array([0, 2, -2, 0, 0, 2, 2, -2, -2])
FM1_vals = np.array([0, 0, 0, 2, -2, 2, -2, 2, -2])
sigma_vals = np.array([10, 50, 50, 70, 70, 100, 40, 40, 100])  # Replace with actual MC results

# Define quadratic function
def cross_section_model(X, c0, c1, c01, sigma_SM):
    FM0, FM1 = X
    return sigma_SM + c0 * FM0**2 + c1 * FM1**2 + c01 * FM0 * FM1

# Initial guess for parameters
initial_guess = [1, 1, 1, 10]

# Fit the function to Monte Carlo data
params, covariance = curve_fit(cross_section_model, (FM0_vals, FM1_vals), sigma_vals, p0=initial_guess)

# Extract fitted coefficients
c0_fit, c1_fit, c01_fit, sigma_SM_fit = params
print(f"Fitted Parameters:\nc0 = {c0_fit}\nc1 = {c1_fit}\nc01 = {c01_fit}\nsigma_SM = {sigma_SM_fit}")
