import numpy as np
from scipy.optimize import root_scalar
import matplotlib.pyplot as plt

# ==========================
# Experiment Setup
# ==========================

luminosity_fb = 100.0                            # [fb^-1]
signal_eff = 20.54 / 100.0                       # Signal efficiency
background_eff = 16.92 / 100.0                   # Background efficiency
sigma_background_fb = 9.9465                     # SM background cross section [fb]
systematic_uncertainty_fraction = 0.10           # 10% systematics

# ==========================
# Derived Quantities
# ==========================

# Number of background events
N_bkg = sigma_background_fb * background_eff * luminosity_fb

# Total uncertainty (statistical + systematic)
delta_stat = np.sqrt(N_bkg)
delta_sys = systematic_uncertainty_fraction * N_bkg
delta_tot = np.sqrt(delta_stat**2 + delta_sys**2)

# Cross-section function in fb for given FM2 (converted from pb)
def sigma_fm2_fb(fm2):
    a = -8.777e-3   # [fb / TeV^4]
    b = 0.5215      # [fb / TeV^8]
    c = 9.941       # [fb]
    return a * fm2 + b * fm2**2 + c

# Number of signal events after selection
def N_sig_fm2(fm2):
    return sigma_fm2_fb(fm2) * signal_eff * luminosity_fb

# Chi^2 definition (background-only expected)
def chi2_fm2(fm2):
    return (N_sig_fm2(fm2) / delta_tot) ** 2

# Solve for FM2 limits (where chi^2 crosses 3.84)
def find_fm2_limits():
    try:
        fm2_upper = root_scalar(lambda f: chi2_fm2(f) - 3.84, bracket=(0.0001, 2.0)).root
        fm2_lower = root_scalar(lambda f: chi2_fm2(f) - 3.84, bracket=(-2.0, -0.0001)).root
        return fm2_lower, fm2_upper
    except Exception:
        return None, None

# Run limit finding
fm2_lower, fm2_upper = find_fm2_limits()

# ==========================
# Print Summary
# ==========================

print("========== 95% CL Exclusion (χ² Method) ==========")
print(f"Background Events (γγ → WW): {N_bkg:.2f}")
print(f"Total Uncertainty (stat + sys): {delta_tot:.2f}")
print(f"Excluded FM2 range at 95% CL: {fm2_lower:.4f} to {fm2_upper:.4f}  [TeV^-4]" if fm2_lower and fm2_upper else "❌ No valid FM2 solution found in range.")

# ==========================
# Plot χ² vs FM2
# ==========================

fm2_vals = np.linspace(-0.05, 0.05, 500)
chi2_vals = [chi2_fm2(f) for f in fm2_vals]

plt.figure(figsize=(8, 5))
plt.plot(fm2_vals, chi2_vals, label=r'$\chi^2(f_{M2})$', color='blue')
plt.axhline(3.84, color='red', linestyle='--', label='95% CL Threshold')
if fm2_lower and fm2_upper:
    plt.axvline(fm2_lower, color='gray', linestyle=':', label=f'FM2 ∈ [{fm2_lower:.4f}, {fm2_upper:.4f}]')
    plt.axvline(fm2_upper, color='gray', linestyle=':')
plt.xlabel(r'$f_{M2}$ [$\mathrm{TeV}^{-4}$]')
plt.ylabel(r'$\chi^2$')
plt.title(r'95% CL Exclusion on $f_{M2}$ (Background-Only)')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()


