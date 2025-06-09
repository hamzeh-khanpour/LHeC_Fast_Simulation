
import numpy as np

# ================================================
# 95% C.L. Exclusion Limit Calculation for FM2
# ================================================

# Experiment Setup
luminosity_fb = 100.0
luminosity_pb = luminosity_fb * 1000.0         # Convert to pb^-1
signal_efficiency = 14.73 / 100.0                 # Final signal efficiency
background_efficiency = 4.46 / 100.0           # Final background efficiency
sigma_background_pb = 0.0099465                # SM WW cross-section [pb]


# Compute expected number of background events
N_b =  sigma_background_pb * luminosity_pb * background_efficiency
print(f"Expected background events (N_b): {N_b:.4f}")

# Compute 95% C.L. upper limit on signal events (Asimov approximation)
N_s_max = 1.96 * np.sqrt(N_b)
print(f"Maximum allowed signal events (N_s^95): {N_s_max:.4f}")

# Convert to cross-section upper limit
sigma_s_max = N_s_max / (luminosity_pb * signal_efficiency)
print(f"Maximum allowed signal cross-section: {sigma_s_max:.6f} pb")

# Quadratic fit parameters (from previous fit.py output)
a = 5.221e-7    # coefficient of FM2^2 (pb/(e-12)^2)
b = 9.942e-3    # SM baseline cross-section (pb)

# Invert the fit to find FM2 value corresponding to sigma_s_max
x_squared = (sigma_s_max - b) / a


if x_squared > 0:
    FM2_limit = np.sqrt(x_squared)
    print(f"95% CL exclusion limit on |FM2|: {FM2_limit:.1f} × 10⁻¹² GeV⁻⁴ = {FM2_limit:.1f} TeV⁻⁴")
else:
    print("❌ No exclusion possible within current fit range: signal cross-section too small.")
