
import numpy as np
from scipy.optimize import fsolve

# ================================================
# 95% C.L. Exclusion and 5σ Discovery for FM2
# Using: σ(FM2) = a * x + b * x² + c
# ================================================

# Inputs
luminosity_fb = 100.0
luminosity_pb = luminosity_fb * 1000.0
signal_efficiency = 14.73 / 100.0
background_efficiency = 4.46 / 100.0
sigma_background_pb = 0.0099465  # pb

# Updated symmetric polynomial fit coefficients
a = -8.615266e-06
b = 5.220918e-07
c = 9.941691e-03

# Background event estimate
N_b = sigma_background_pb * luminosity_pb * background_efficiency
N_s_95 = 1.96 * np.sqrt(N_b)
sigma_s_max = N_s_95 / (luminosity_pb * signal_efficiency)



# Solve: a x + b x^2 + c = sigma_s_max
A = b
B = a
C = c - sigma_s_max

discriminant = B**2 - 4 * A * C

print(f"Expected background events (N_b): {N_b:.4f}")
print(f"Maximum allowed signal events (N_s^95): {N_s_95:.4f}")
print(f"Maximum allowed signal cross-section: {sigma_s_max:.6f} pb")

if discriminant >= 0 and A != 0:
    x1 = (-B + np.sqrt(discriminant)) / (2 * A)
    x2 = (-B - np.sqrt(discriminant)) / (2 * A)
    roots = [r for r in [x1, x2] if r > 0]
    FM2_excl = min(roots) if roots else None
    if FM2_excl:
        print(f"✅ 95% CL exclusion limit on |FM2|: {FM2_excl:.1f} TeV⁻⁴")
    else:
        print("❌ No positive solution for FM2 exclusion.")
else:
    print("❌ No exclusion possible within current fit range.")




# === Discovery: 5σ threshold ===
def solve_discovery_s(N_b, significance=5.0):
    f = lambda s: s / np.sqrt(s + N_b) - significance
    return fsolve(f, significance * np.sqrt(N_b))[0]

N_s_disc = solve_discovery_s(N_b)
sigma_s_disc = N_s_disc / (luminosity_pb * signal_efficiency)

# Solve: a x + b x^2 + c = sigma_s_disc
C_disc = c - sigma_s_disc
disc_discriminant = B**2 - 4 * A * C_disc

print(f"\n🔍 5σ Discovery Threshold:")
if disc_discriminant >= 0 and A != 0:
    x1 = (-B + np.sqrt(disc_discriminant)) / (2 * A)
    x2 = (-B - np.sqrt(disc_discriminant)) / (2 * A)
    roots_disc = [r for r in [x1, x2] if r > 0]
    FM2_disc = min(roots_disc) if roots_disc else None
    if FM2_disc:
        print(f"  |FM2| > {FM2_disc:.1f} TeV⁻⁴")
    else:
        print("❌ No valid discovery solution — signal too weak")
else:
    print("❌ No real solution — discovery cross-section too small")
