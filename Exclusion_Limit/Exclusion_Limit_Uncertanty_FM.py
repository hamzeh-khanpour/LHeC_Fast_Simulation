# Hamzeh Khanpour

import numpy as np
import scipy.stats as stats
import scipy.optimize as opt
import scipy.integrate as integrate
import math
import matplotlib.pyplot as plt
from scipy.stats import poisson
import mplhep as hep

# Matplotlib configuration for publication-quality plots
hep.style.use("CMS")


# Cross-section values (in pb)
sigma_SM    = 0.0134936  # Standard Model cross-section

sigma_FM0_1 = 1.99216    # Cross-section for FM0 = 10^-8
sigma_FM1_1 = 0.17783    # Cross-section for FM1 = 10^-8
sigma_FM2_1 = 77.8809    # Cross-section for FM2 = 10^-8
sigma_FM3_1 = 5.95386    # Cross-section for FM3 = 10^-8


# Given parameters
luminosity = 0.10  # in pb^-1


# FM0 value from anoinputs file in MG
FM0_value = 10**-8  # The input value for FM0 in the cross-section calculation in MG.


# Compute A coefficient for anomalous contribution with correct scaling
A = (sigma_FM0_1 - sigma_SM) / FM0_value**2  # Proper scaling for FM0


signal_efficiency = 26.11 / 100.0  # Efficiency for signal as decimal
background_efficiency = 0.99 / 100.0  # Efficiency for background as decimal


# Compute expected background count
background_events = background_efficiency * sigma_SM * luminosity


# Systematic uncertainty (assumed 10%)
delta_sys = 0.10


# Compute global error term (must be outside function)
N_SM = luminosity * sigma_SM
delta_st = 1 / np.sqrt(max(N_SM, 1e-10))  # Statistical uncertainty with floor value
error_term = sigma_SM * np.sqrt(delta_st**2 + delta_sys**2)  # Total uncertainty


# Poisson probability function
def poisson_prob(n, S, B):
    """
    Poisson probability of observing n events given S (signal) and B (background).
    Uses scipy.stats.poisson for numerical stability.
    """
    mean = S + B
    return poisson.pmf(n, mean)


# Compute Bayesian upper limit for sigma_CL using integration
def compute_sigma_CL(CL, esig, L, B):
    """
    Bayesian confidence level calculation for the upper bound on signal cross-section,
    including both statistical and systematic uncertainties.
    """
    n = round(B)

    def integrand(osig):
        return poisson_prob(n, osig * esig * L + error_term, B)  # Now includes uncertainty

    # Define integration range for normalization
    upper_bound = max(100 * sigma_SM, 500 * sigma_SM, 1000 * sigma_SM)

    den, _ = integrate.quad(integrand, 0, upper_bound)

    def find_sigma_CLi(sigma_CLi):
        num, _ = integrate.quad(integrand, 0, sigma_CLi)
        return num - CL * den

    try:
        result = opt.root_scalar(find_sigma_CLi, bracket=[0, upper_bound], method='brentq')
        return result.root if result.converged else None
    except Exception as e:
        print(f"Error solving for sigma_CL: {e}")
        return None

# Compute 95% CL upper limit on sigma_S
sigma_CL = compute_sigma_CL(0.95, signal_efficiency, luminosity, background_events)


# Compute F_i exclusion limit using the Bayesian upper limit on sigma_S
if sigma_CL > sigma_SM + error_term and (sigma_CL - sigma_SM - error_term) > 0:
    F_i_limit = np.sqrt((sigma_CL - sigma_SM - error_term) / A)
else:
    F_i_limit = 0  # Prevent NaN or negative sqrt


# Print Final Results
print("\n========= FINAL RESULTS =========")
print(f"Computed A coefficient: {A:.4e}")
print(f"📌 95% CL upper limit on cross-section: {sigma_CL:.6f} pb")
print(f"📌 95% CL exclusion limit on F_M0: {F_i_limit:.4e}")
print("==================================\n")


# Compute expected signal events for F_i values
if F_i_limit > 0:
    F_i_values = np.linspace(0, 2 * F_i_limit, 1000)
else:
    F_i_values = np.linspace(0, 1e-7, 1000)  # Small range to avoid flat plots

N_S_values = (sigma_SM + A * F_i_values**2) * luminosity * signal_efficiency


plt.figure(figsize=(9, 10))
plt.subplots_adjust(left=0.15, right=0.95, bottom=0.12, top=0.95)

plt.plot(F_i_values, N_S_values, label=r'$N_S(f_{M_0} / \Lambda^4)$', linewidth=3)
plt.axhline(sigma_CL * luminosity * signal_efficiency, color='r', linestyle='--', label=r'$S_{95}$ limit', linewidth=3)
plt.axvline(F_i_limit, color='g', linestyle='--', label=r'$(f_{M_0} / \Lambda^4)^{lim}$', linewidth=3)

plt.xlabel(r'95% CL on $f_{M_0} / \Lambda^4$', labelpad=15, loc='center')
plt.ylabel(r'Expected Signal Events $N_S$', labelpad=15, loc='center')
plt.legend()
plt.title(r'Preliminary : LHeC@1.2 TeV (100 fb$^{-1}$)')
plt.grid()

# Save the figure as a JPG file  
plt.savefig("LHeC_Limit_FM0.jpg", dpi=300, bbox_inches='tight', format='jpg')

plt.show()
