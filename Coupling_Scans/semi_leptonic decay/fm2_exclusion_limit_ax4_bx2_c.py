import numpy as np

# ================================================
# 95% C.L. Exclusion Limit Calculation for FM2
# Using Fit: σ = a x⁴ + b x² + c
# ================================================

# Experiment Setup
luminosity_fb = 100.0
luminosity_pb = luminosity_fb * 1000.0         # Convert to pb^-1
signal_efficiency = 2.33 / 100.0                     # Final signal efficiency
background_efficiency = 0.41 / 100.0                  # Final background efficiency
sigma_background_pb = 0.0099465                # SM WW cross-section [pb]

# Compute expected number of background events
N_b = sigma_background_pb * luminosity_pb * background_efficiency
print(f"Expected background events (N_b): {N_b:.4f}")

# Compute 95% CL upper limit on signal events (Asimov approximation)
N_s_max = 1.96 * np.sqrt(N_b)
print(f"Maximum allowed signal events (N_s^95): {N_s_max:.4f}")

# Convert to cross-section upper limit
sigma_s_max = N_s_max / (luminosity_pb * signal_efficiency)
print(f"Maximum allowed signal cross-section: {sigma_s_max:.6f} pb")


# Quartic + quadratic fit parameters: σ = a x⁴ + b x² + c
a = -1.361928e-17  # coefficient of FM2^4
b = 5.221039e-7    # coefficient of FM2^2
c = 9.941209e-3    # SM cross-section [pb]

# Solve a x^4 + b x^2 + c = sigma_s_max → a y^2 + b y + (c - sigma_s_max) = 0, with y = x^2
A = a
B = b
C = c - sigma_s_max

discriminant = B**2 - 4 * A * C

if A != 0 and discriminant >= 0:
    y1 = (-B + np.sqrt(discriminant)) / (2 * A)
    y2 = (-B - np.sqrt(discriminant)) / (2 * A)

    FM2_limit_1 = np.sqrt(y1) if y1 >= 0 else None
    FM2_limit_2 = np.sqrt(y2) if y2 >= 0 else None

    print("\n95% CL exclusion limit on |FM2| (real positive roots only):")
    if FM2_limit_1:
        print(f"  |FM2| < {FM2_limit_1:.1f} × 10⁻¹² GeV⁻⁴ = {FM2_limit_1:.1f} TeV⁻⁴")
    if FM2_limit_2 and (not FM2_limit_1 or abs(FM2_limit_2) < abs(FM2_limit_1)):
        print(f"  |FM2| < {FM2_limit_2:.1f} × 10⁻¹² GeV⁻⁴ = {FM2_limit_2:.1f} TeV⁻⁴")
    if not FM2_limit_1 and not FM2_limit_2:
        print("❌ No real, positive solutions for FM2.")
else:
    print("❌ No exclusion possible within current fit range: signal cross-section too small.")

