
import numpy as np

# === Setup ===
luminosity_fb = 100.0  # in fb⁻¹
luminosity_pb = luminosity_fb * 1000.0  # in fb⁻¹     ?
signal_efficiency = 14.73 / 100.0
background_efficiency = 4.46 / 100.0
sigma_background_pb = 0.0099465  # pb

# Symmetric polynomial fit: σ(FM2) = a*x + b*x^2 + c
a = -8.615266e-06
b = 5.220918e-07
c = 9.941691e-03

# Compute expected number of background events
N_b = sigma_background_pb * luminosity_pb * background_efficiency

# Compute 95% C.L. upper limit on signal events
N_s_max = 1.96 * np.sqrt(N_b)

# Convert to signal cross-section upper limit
sigma_s_max = N_s_max / (luminosity_pb * signal_efficiency)

# Solve: a x + b x^2 + c = sigma_s_max --> b x^2 + a x + (c - sigma_s_max) = 0
A = b
B = a
C = c - sigma_s_max

discriminant = B**2 - 4 * A * C
FM2_excl = None
if discriminant >= 0 and A != 0:
    x1 = (-B + np.sqrt(discriminant)) / (2 * A)
    x2 = (-B - np.sqrt(discriminant)) / (2 * A)
    roots = [r for r in [x1, x2] if r > 0]
    FM2_excl = min(roots) if roots else None

# Output result
print(f"Expected Background Events (N_b): {N_b:.4f}")
print(f"Maximum allowed signal cross-section: {sigma_s_max:.6f} pb")
if FM2_excl:
    print(f"✅ 95% C.L. exclusion limit on |FM2|: {FM2_excl:.1f} TeV⁻⁴")
else:
    print("❌ No exclusion possible — signal cross-section too small")
