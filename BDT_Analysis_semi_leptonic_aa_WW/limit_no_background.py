# ================================================================================
#   Hamzeh Khanpour — June 2025
#   No-Background 95% CL Limit on FM2 (cut-and-count + profile likelihood)
# ================================================================================

import numpy as np
from scipy.optimize import root_scalar
import matplotlib.pyplot as plt

# --- Inputs ---
luminosity_fb = 100.0
delta_sys = 0.10  # 10% systematic uncertainty

signal_efficiency = 0.1405  # Example from your latest BDT output
background_efficiency  = 0.1305

sigma_sm_fb = 1.538327e+01 * background_efficiency # SM cross section in fb

# --- Cross-section function (in fb) ---
def sigma_fm2_fb(fm2):
    a = -1.464725e-02
    b = 8.929319e-03
    c = 1.538327e+01
    return a * fm2 + b * fm2**2 + c

# --- Derived ---
n_sm = luminosity_fb * sigma_sm_fb
delta_stat = 1.0 / np.sqrt(n_sm)
delta_tot = np.sqrt(delta_stat**2 + delta_sys**2)

# --- χ² test function ---
def chi2_stat(fm2):
    sigma_bsm = sigma_fm2_fb(fm2) * signal_efficiency
    return ((sigma_sm_fb - sigma_bsm) / (sigma_sm_fb * delta_tot)) ** 2

# --- Find limits where χ² = 3.84 ---
def chi2_diff(fm2):
    return chi2_stat(fm2) - 3.84

# Scan around fm2 = 0
result_lower = root_scalar(chi2_diff, bracket=[-50, 0.0001], method='brentq')
result_upper = root_scalar(chi2_diff, bracket=[0.0001, 50], method='brentq')

# --- Output ---
print("========= χ² Limit Setting at 95% CL =========")
print(f"Luminosity: {luminosity_fb} fb⁻¹")
print(f"SM cross section: {sigma_sm_fb:.4f} fb")
print(f"Stat. error: {delta_stat:.4f}, Total rel. error: {delta_tot:.4f}")
print(f"95% CL Exclusion Range for fM2: {result_lower.root:.3f} to {result_upper.root:.3f} [TeV⁻⁴]")

# --- Optional: Plot χ² curve ---
fm2_vals = np.linspace(-50, 50, 200)
chi2_vals = [chi2_stat(fm2) for fm2 in fm2_vals]

plt.plot(fm2_vals, chi2_vals, label=r"$\chi^2(f_{M2})$")
plt.axhline(3.84, color='red', linestyle='--', label='95% CL Threshold')
plt.xlabel(r"$f_{M2} / \Lambda^4\ [\mathrm{TeV}^{-4}]$")
plt.ylabel(r"$\chi^2$")
plt.legend()
plt.grid(True)
plt.title(r"Profile $\chi^2$ for $f_{M2}$ at 95% CL")
plt.show()
