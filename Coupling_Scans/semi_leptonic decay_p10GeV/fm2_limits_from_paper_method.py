
import numpy as np
from scipy.optimize import fsolve
import matplotlib.pyplot as plt

# === Physics Parameters from Paper (Asimov-style method) ===
luminosity_fb = 100.0
luminosity_pb = luminosity_fb * 1000.0
signal_eff = 14.73 / 100.0
background_eff = 4.46 / 100.0
sigma_sm_pb = 0.0099465  # SM cross-section in pb

# Fit model: σ(fM2) = a * fM2^2 + b
a = 5.221e-7  # pb / (TeV^-4)^2
b = sigma_sm_pb  # SM cross-section

# Expected background
N_b = sigma_sm_pb * background_eff * luminosity_pb

# === 1. 95% C.L. Exclusion Limit using Asimov approximation ===
N_s_95 = 1.96 * np.sqrt(N_b)
sigma_95 = N_s_95 / (luminosity_pb * signal_eff)

# Invert fit: σ = a * x^2 + b → x^2 = (σ - b)/a
x2_excl = (sigma_95 - b) / a
fM2_excl = np.sqrt(x2_excl) if x2_excl > 0 else None

# === 2. 5σ Discovery Limit using s / sqrt(s + b) = 5 ===
def discovery_significance_equation(s, B, z=5.0):
    return s / np.sqrt(s + B) - z

N_s_disc = fsolve(discovery_significance_equation, x0=5*np.sqrt(N_b), args=(N_b))[0]
sigma_disc = N_s_disc / (luminosity_pb * signal_eff)
x2_disc = (sigma_disc - b) / a
fM2_disc = np.sqrt(x2_disc) if x2_disc > 0 else None

# === Summary Output ===
print(f"Expected Background Events (N_b): {N_b:.4f}")
print(f"95% CL Signal Events (N_s^95): {N_s_95:.4f}")
print(f"95% CL Max σ_s: {sigma_95:.6f} pb")
print(f"5σ Discovery σ_s: {sigma_disc:.6f} pb")
print()
if fM2_excl:
    print(f"✅ 95% CL Exclusion Limit: |FM2| < {fM2_excl:.1f} TeV⁻⁴")
else:
    print("❌ No exclusion possible — cross-section too small.")
if fM2_disc:
    print(f"✅ 5σ Discovery Threshold: |FM2| > {fM2_disc:.1f} TeV⁻⁴")
else:
    print("❌ No discovery threshold — signal too small.")

# === Plotting ===
fM2_range = np.linspace(0, 1000, 500)
sigma_fM2 = a * fM2_range**2 + b

plt.figure(figsize=(8, 5))
plt.plot(fM2_range, sigma_fM2, label=r"$\sigma(FM2) = a x^2 + b$")
plt.axhline(sigma_95, color='red', linestyle='--', label=r"95% CL Exclusion $\sigma$")
plt.axhline(sigma_disc, color='orange', linestyle='--', label=r"5$\sigma$ Discovery $\sigma$")
plt.xlabel(r"$|FM2|$ [$\mathrm{TeV}^{-4}$]")
plt.ylabel(r"Cross-section [pb]")
plt.title("Cross-Section vs |FM2| with 95% Exclusion and 5σ Discovery")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig("fm2_limits_from_paper_method.pdf")
plt.show()
