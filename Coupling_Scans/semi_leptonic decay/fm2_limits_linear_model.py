
import numpy as np
from scipy.optimize import fsolve

# ================================================
# 95% C.L. Exclusion Limit Calculation for FM2 (Linear σ model)
# ================================================

# Experiment Setup
luminosity_fb = 100.0
luminosity_pb = luminosity_fb * 1000.0         # Convert to pb^-1
signal_efficiency = 14.73 / 100.0              # Final signal efficiency
background_efficiency = 4.46 / 100.0           # Final background efficiency
sigma_sm_pb = 0.0099465                        # SM WW cross-section [pb]
A = 5.221e-7                                   # linear coefficient: sigma = sigma_SM + A * FM2

# Step 1: Compute expected background yield
N_b = sigma_sm_pb * background_efficiency * luminosity_pb
print(f"Expected background events (N_b): {N_b:.4f}")

# Step 2: Asimov 95% CL exclusion
# Solve: effS * (sigma_SM + A * FM2) * L / sqrt(N_b) = 1.95
numerator = 1.95 * np.sqrt(N_b)
denominator = signal_efficiency * A * luminosity_pb
FM2_excl = (numerator / denominator) - (sigma_sm_pb / A)

print(f"✅ 95% CL exclusion limit on |FM2|: {FM2_excl:.1f} TeV⁻⁴")

# Step 3: 5σ Discovery Threshold
def discovery_significance_equation(s, B, z=5.0):
    return s / np.sqrt(s + B) - z

N_s_disc = fsolve(discovery_significance_equation, x0=5*np.sqrt(N_b), args=(N_b))[0]
sigma_disc = N_s_disc / (luminosity_pb * signal_efficiency)
FM2_disc = (sigma_disc - sigma_sm_pb) / A

print("🔍 5σ Discovery Threshold:")
if FM2_disc > 0:
    print(f"  |FM2| > {FM2_disc:.1f} TeV⁻⁴")
else:
    print("❌ No valid solution — signal too small")
