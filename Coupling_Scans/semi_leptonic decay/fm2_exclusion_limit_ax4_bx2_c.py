import numpy as np
from scipy.optimize import fsolve

# ================================================
# 95% C.L. Exclusion and 5σ Discovery for FM2
# ================================================

# Inputs
luminosity_fb = 100.0
luminosity_pb = luminosity_fb * 1000.0
signal_efficiency = 14.73 / 100.0
background_efficiency = 4.46 / 100.0
sigma_background_pb = 0.0099465  # pb

# Quartic fit coefficients: σ(FM2) = a x⁴ + b x² + c
a = -1.361928e-17
b = 5.221039e-7
c = 9.941209e-3

# === Exclusion: compute N_b and limit ===
N_b = sigma_background_pb * luminosity_pb * background_efficiency
N_s_95 = 1.96 * np.sqrt(N_b)
sigma_s_max = N_s_95 / (luminosity_pb * signal_efficiency)

# Invert quartic: let y = x², solve a y² + b y + (c - sigma_s_max) = 0
A = a
B = b
C = c - sigma_s_max
discriminant = B**2 - 4 * A * C

print(f"Expected background events (N_b): {N_b:.4f}")
print(f"Maximum allowed signal events (N_s^95): {N_s_95:.4f}")
print(f"Maximum allowed signal cross-section: {sigma_s_max:.6f} pb")

if discriminant >= 0 and A != 0:
    y1 = (-B + np.sqrt(discriminant)) / (2 * A)
    y2 = (-B - np.sqrt(discriminant)) / (2 * A)

    FM2_excl_1 = np.sqrt(y1) if y1 > 0 else None
    FM2_excl_2 = np.sqrt(y2) if y2 > 0 else None

    print("\n✅ 95% CL exclusion limit on |FM2|:")
    if FM2_excl_1:
        print(f"  |FM2| < {FM2_excl_1:.1f} TeV⁻⁴")
    if FM2_excl_2 and (not FM2_excl_1 or FM2_excl_2 < FM2_excl_1):
        print(f"  |FM2| < {FM2_excl_2:.1f} TeV⁻⁴")
    if not FM2_excl_1 and not FM2_excl_2:
        print("❌ No real positive solution for FM2.")
else:
    print("❌ No exclusion possible within current fit range.")

# === Discovery: solve for N_s such that s / sqrt(s + b) = 5 ===
def solve_discovery_s(N_b, significance=5.0):
    f = lambda s: s / np.sqrt(s + N_b) - significance
    return fsolve(f, significance * np.sqrt(N_b))[0]

N_s_disc = solve_discovery_s(N_b)
sigma_s_disc = N_s_disc / (luminosity_pb * signal_efficiency)

# Invert quartic fit for discovery cross-section
C_disc = c - sigma_s_disc
discriminant_disc = B**2 - 4 * A * C_disc

print(f"\n🔍 5σ Discovery Threshold:")
if discriminant_disc >= 0 and A != 0:
    y1 = (-B + np.sqrt(discriminant_disc)) / (2 * A)
    y2 = (-B - np.sqrt(discriminant_disc)) / (2 * A)

    FM2_disc_1 = np.sqrt(y1) if y1 > 0 else None
    FM2_disc_2 = np.sqrt(y2) if y2 > 0 else None

    if FM2_disc_1:
        print(f"  |FM2| > {FM2_disc_1:.1f} TeV⁻⁴")
    if FM2_disc_2 and (not FM2_disc_1 or FM2_disc_2 < FM2_disc_1):
        print(f"  |FM2| > {FM2_disc_2:.1f} TeV⁻⁴")
    if not FM2_disc_1 and not FM2_disc_2:
        print("❌ No valid discovery solution — signal too weak")
else:
    print("❌ No real solution — discovery cross-section too small")
