import numpy as np
import scipy.stats as stats
import scipy.optimize as opt
import scipy.integrate as integrate
import math
import matplotlib.pyplot as plt
from scipy.stats import poisson

# Given parameters
luminosity = 10  # in ab^-1
signal_efficiency = 0.5  # Efficiency of detecting signal
background_events = 100  # Expected background count

# Cross-section values (in pb)
sigma_SM = 0.0134936  # Standard Model cross-section
sigma_FM0_1 = 1.99216  # Cross-section for FM0 = 1

# Compute A coefficient for anomalous contribution
A = sigma_FM0_1 - sigma_SM

# Systematic uncertainty (assumed 5%)
delta_sys = 0.05

# Use Scipy's poisson function to avoid factorial overflow
def poisson_prob(n, S, B):
    """
    Poisson probability of observing n events given S (signal) and B (background).
    Uses scipy.stats.poisson for numerical stability.
    """
    mean = S + B
    return poisson.pmf(n, mean)  # More stable calculation

# Compute Bayesian upper limit for sigma_CL using integration
def compute_sigma_CL(CL, esig, L, B):
    """
    Bayesian confidence level calculation for the upper bound on signal cross-section.
    Equivalent to Mathematica's σCL function, now including statistical and systematic uncertainties.
    """
    n = round(B)  # Convert background count to integer
    N_SM = L * sigma_SM  # Compute expected SM events
    delta_st = 1 / np.sqrt(N_SM)  # Statistical uncertainty
    error_term = sigma_SM * np.sqrt(delta_st**2 + delta_sys**2)  # Total uncertainty

    # Compute denominator (full probability normalization)
    def integrand(osig):
        return poisson_prob(n, osig * esig * L, B)

    upper_bound = 5 * sigma_FM0_1  # Avoid infinite integration issues
    den, _ = integrate.quad(integrand, 0, upper_bound)  # Fixed integration limit

    # Find sigma_CL by integrating up to σCLi
    def find_sigma_CLi(sigma_CLi):
        num, _ = integrate.quad(integrand, 0, sigma_CLi)  # Fixed integration
        return num - CL * den  # Solve num = CL * den

    # Solve for sigma_CL numerically
    result = opt.root_scalar(find_sigma_CLi, bracket=[0, upper_bound], method='brentq')
    return result.root if result.converged else None

# Compute 95% CL upper limit on sigma_S
sigma_CL = compute_sigma_CL(0.95, signal_efficiency, luminosity, background_events)

# Compute F_i exclusion limit using the Bayesian upper limit on sigma_S
F_i_limit = np.sqrt((sigma_CL - sigma_SM) / A) if sigma_CL > sigma_SM else 0

print(f"95% CL Bayesian upper limit on cross-section (with uncertainties): {sigma_CL:.6f} pb")
print(f"95% CL exclusion limit on F_i: {F_i_limit:.4f}")

# Plotting
F_i_values = np.linspace(0, 2 * F_i_limit, 100)
N_S_values = (sigma_SM + A * F_i_values**2) * luminosity * signal_efficiency

plt.figure(figsize=(8, 6))
plt.plot(F_i_values, N_S_values, label=r'$N_S(F_i)$')
plt.axhline(sigma_CL * luminosity * signal_efficiency, color='r', linestyle='--', label=r'$S_{95}$ limit')
plt.axvline(F_i_limit, color='g', linestyle='--', label=r'$F_i^{lim}$')
plt.xlabel(r'$F_i$')
plt.ylabel(r'Expected Signal Events $N_S$')
#plt.yscale('log')  # Log scale for better visualization
plt.legend()
plt.title('95% CL Bayesian Limit on $F_i$')
plt.grid()
plt.show()
