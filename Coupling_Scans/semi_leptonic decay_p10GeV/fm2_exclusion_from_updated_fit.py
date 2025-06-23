
import numpy as np

# ================================================
# 95% C.L. Exclusion Limit Calculation for FM2
# Using: σ(FM2) = a * FM2 + b * FM2² + c
# ================================================

# Inputs
luminosity_fb = 100.0
luminosity_pb = luminosity_fb * 1000.0
signal_eff = 14.73 / 100.0
background_eff = 4.46 / 100.0
sigma_bkg = 0.0099465  # pb


# Fit coefficients
a = -8.777e-06
b = 5.215e-07
c = 9.941e-03


# Step 1: Compute expected background events
N_b = sigma_bkg * background_eff * luminosity_pb

# Step 2: Compute 95% CL upper limit on signal events
N_s_max = 1.96 * np.sqrt(N_b)

# Step 3: Convert to cross-section
sigma_s_max = N_s_max / (luminosity_pb * signal_eff)

# Step 4: Solve a x + b x^2 + c = sigma_s_max => b x^2 + a x + (c - sigma_s_max) = 0
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

# Print results
print(f"Expected Background Events (N_b): {N_b:.4f}")
print(f"Max allowed signal cross-section (σ_s^max): {sigma_s_max:.6f} pb")
if FM2_excl:
    print(f"✅ 95% CL exclusion limit on |FM2|: {FM2_excl:.1f} TeV⁻⁴")
else:
    print("❌ No exclusion possible — signal cross-section too small.")
