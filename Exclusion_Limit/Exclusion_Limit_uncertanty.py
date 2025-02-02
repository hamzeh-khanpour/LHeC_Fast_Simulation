
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


# Cross-section values (in pb -> converted to fb)
sigma_SM = 0.0134936 * 1000   # Standard Model cross-section
sigma_FM0_1 = 1.99385 * 1000  # Cross-section for FM0 = 10^-8

# Given parameters
luminosity = 100.0  # in fb^-1

signal_efficiency = 49.11 / 100.0  # Efficiency as decimal
background_efficiency = 0.99 / 100.0  # Efficiency as decimal

# Compute expected background count
background_events = background_efficiency * sigma_SM * luminosity

# Systematic uncertainty (assumed 5%)
delta_sys = 0.05

# FM0 value from anoinputs file
FM0_value = 1.0 # 10**-8  # Given in anoinputs file

# Compute A coefficient for anomalous contribution with correct scaling
A = (sigma_FM0_1 - sigma_SM) / FM0_value**2  # Proper scaling for FM0

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
    Bayesian confidence level calculation for the upper bound on signal cross-section.
    """
    n = round(B)
    N_SM = L * sigma_SM
    delta_st = 1 / np.sqrt(N_SM)  # Statistical uncertainty
    error_term = sigma_SM * np.sqrt(delta_st**2 + delta_sys**2)  # Total uncertainty

    def integrand(osig):
        return poisson_prob(n, osig * esig * L, B)

    upper_bound = max(10 * sigma_SM, 5 * sigma_FM0_1)
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
F_i_limit = np.sqrt((sigma_CL - sigma_SM) / A) if sigma_CL > sigma_SM else 0

print(f"95% CL Bayesian upper limit on cross-section (with uncertainties): {sigma_CL:.6f} pb")
print(f"95% CL exclusion limit on F_i: {F_i_limit:.4f}")

# Plotting
F_i_values = np.linspace(0, 2 * F_i_limit, 100)
N_S_values = (sigma_SM + A * F_i_values**2) * luminosity * signal_efficiency

plt.subplots_adjust(left=0.15, right=0.95, bottom=0.12, top=0.95)
plt.plot(F_i_values, N_S_values, label=r'$N_S(F_i)$')
plt.axhline(sigma_CL * luminosity * signal_efficiency, color='r', linestyle='--', label=r'$S_{95}$ limit')
plt.axvline(F_i_limit, color='g', linestyle='--', label=r'$F_i^{lim}$')
plt.xlabel(r'$F_i$')
plt.ylabel(r'Expected Signal Events $N_S$')
plt.legend()
plt.title('95% CL Bayesian Limit on $F_i$')
plt.grid()
plt.show()


