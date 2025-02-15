import numpy as np

# Given values from table
sigma_SM = 0.0134936  # pb
luminosity = 100  # fb^-1
eff_background = 0.0099  # 0.99% background efficiency (decimal)
eff_signal = np.mean([0.2632, 0.2332, 0.2607, 0.2603, 0.2672, 0.2599, 0.2664, 0.2646])  # Average signal efficiency

# Compute expected background events
N_bkg = sigma_SM * luminosity * eff_background

# Compute uncertainties
delta_stat = np.sqrt(N_bkg)
delta_sys = 0.10 * N_bkg
delta_total = np.sqrt(delta_stat**2 + delta_sys**2)

# Compute upper limit on signal events at 95% CL
N_sig_limit = N_bkg + 1.96 * delta_total

# Compute sigma_limit
sigma_limit = N_sig_limit / (luminosity * eff_signal/100.0)  # Convert to pb

print(f"95% CL Upper Limit on Cross-Section: {sigma_limit:.6f} pb")
