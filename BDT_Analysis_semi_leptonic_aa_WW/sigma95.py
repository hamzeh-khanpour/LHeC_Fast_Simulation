import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import poisson

# -----------------------------------
# Inputs (your values)
# -----------------------------------
luminosity_fb = 100.0
luminosity_pb = luminosity_fb * 1000.0  # Convert to pb^-1

sigma_SM = 0.0099465       # pb
sigma_FM2_100 = 0.0142882  # pb for fM2 = 100 TeV^-4

signal_eff = 0.0233        # 2.33%
bkg_eff = 0.0041           # 0.41%

fM2_ref = 100              # Reference value used in simulation

# -----------------------------------
# Calculate σ_BSM
# -----------------------------------
sigma_BSM = (sigma_FM2_100 - sigma_SM) / (fM2_ref ** 2)

# Background event yield
N_bkg = sigma_SM * luminosity_pb * bkg_eff
N_obs = int(np.round(N_bkg))  # Assume observed = expected_bkg
N_95 = poisson.ppf(0.95, N_obs)

# -----------------------------------
# Scan over fM2 in TeV⁻⁴
# -----------------------------------
fM2_vals = np.linspace(-300, 300, 1000)
sigma_fM2 = sigma_SM + (fM2_vals ** 2) * sigma_BSM
N_sig = sigma_fM2 * luminosity_pb * signal_eff
N_tot = N_sig + N_bkg

# Exclusion region
excluded_mask = N_tot > N_95
excluded_fM2_vals = fM2_vals[excluded_mask]

# -----------------------------------
# Plot
# -----------------------------------
plt.figure(figsize=(10, 6))
plt.plot(fM2_vals, N_tot, label='Signal + Background Yield')
plt.axhline(N_95, color='r', linestyle='--', label=f'95% CL Threshold = {N_95:.1f} events')
plt.fill_between(fM2_vals, 0, N_tot, where=excluded_mask, color='orange', alpha=0.3, label='Excluded')
plt.xlabel(r'$f_{M2}$ [$\mathrm{TeV}^{-4}$]')
plt.ylabel('Expected Event Yield')
plt.title('95% CL Exclusion Limit on $f_{M2}$')
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()

# -----------------------------------
# Result
# -----------------------------------
if excluded_fM2_vals.size > 0:
    fM2_excluded_min = excluded_fM2_vals.min()
    fM2_excluded_max = excluded_fM2_vals.max()
    print(f"Excluded region at 95% CL: fM2 ∈ [{fM2_excluded_min:.1f}, {fM2_excluded_max:.1f}] TeV⁻⁴")
else:
    print("No exclusion found at 95% CL.")
