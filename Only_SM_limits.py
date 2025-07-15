import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import root_scalar

# Input parameters
sigma_SM_fb = 15.4074  # Standard Model cross section [fb]
luminosity_fb = 1000.0  # Integrated luminosity [fb^-1]
delta_sys = 0.0   # Systematic uncertainty (10%)

# Calculate statistical uncertainty
N_SM = sigma_SM_fb * luminosity_fb
delta_stat = 1 / np.sqrt(N_SM)

# Define sigma(fM2) function in fb
def sigma_fm2_fb(fm2):
    a = -1.464725e-02
    b = 8.929319e-04
    c = 1.538327e+01
    return a * fm2 + b * fm2**2 + c

# Chi-squared function
def chi2(fm2):
    sigma_BSM = sigma_fm2_fb(fm2)
    numerator = sigma_SM_fb - sigma_BSM
    denominator = sigma_SM_fb * np.sqrt(delta_stat**2 + delta_sys**2)
    return (numerator / denominator) ** 2

# Function to solve for 95% CL limits (χ² = 3.84)
def chi2_minus_3p84(fm2):
    return chi2(fm2) - 3.84

# Find limits numerically
lower_limit = root_scalar(chi2_minus_3p84, bracket=[-50, 0.0001]).root
upper_limit = root_scalar(chi2_minus_3p84, bracket=[0.0001, 50]).root

# Print results
print(f"95% CL limits on fM2/Λ⁴: [{lower_limit:.2f}, {upper_limit:.2f}] TeV⁻⁴")

# Plot cross section vs fM2
fm2_values = np.linspace(-50, 50, 1000)
sigma_values = sigma_fm2_fb(fm2_values)

plt.figure()
plt.plot(fm2_values, sigma_values, label=r'$\sigma(f_{M2})$ [fb]')
plt.axhline(sigma_SM_fb, color='gray', linestyle='--', label=r'$\sigma_{\mathrm{SM}}$')
plt.xlabel(r'$f_{M2} \, [\mathrm{TeV}^{-4}]$')
plt.ylabel(r'$\sigma$ [fb]')
plt.title('Cross Section vs $f_{M2}$')
plt.legend()
plt.grid(True)
plt.show()

# Plot chi-squared vs fM2
chi2_values = chi2(fm2_values)

plt.figure()
plt.plot(fm2_values, chi2_values, label=r'$\chi^2(f_{M2})$')
plt.axhline(3.84, color='red', linestyle='--', label='95% CL')
plt.xlabel(r'$f_{M2} \, [\mathrm{TeV}^{-4}]$')
plt.ylabel(r'$\chi^2$')
plt.title('Chi-Squared vs $f_{M2}$')
plt.legend()
plt.grid(True)
plt.show()
