
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve, brentq
from scipy.stats import poisson

# === Setup ===
luminosity_fb = 100.0
luminosity_pb = luminosity_fb * 1000.0
signal_eff = 14.73 / 100.0
background_eff = 4.46 / 100.0
sigma_sm_pb = 0.0099465

# Fit: σ(fM2) = a * fM2^2 + b
a = 5.221e-7  # pb / (TeV^-4)^2
b = sigma_sm_pb


# Expected background events
N_b = sigma_sm_pb * background_eff * luminosity_pb

# === 1. Bayesian Poisson Limit ===
# Goal: find sigma_s^max such that ∫₀^{N_obs} P(n | s + b) dn = 0.95

# Bayesian one-sided upper limit for Poisson process
def bayesian_limit(N_obs, N_b, CL=0.95):
    def integrand(s, n_obs, bkg):
        lam = s + bkg
        return poisson.pmf(n_obs, lam)

    def cumulative_probability(s_max):
        return sum(poisson.pmf(N_obs, s_max + N_b) for N_obs in range(int(N_b) + 1))

    def root_func(s_max):
        return cumulative_probability(s_max) - CL

    s_max = brentq(root_func, 0, 50)
    return s_max


# Use expected number of events as "observation" in the Asimov sense
N_obs = int(round(N_b))
N_s_bayes = bayesian_limit(N_obs, N_b, CL=0.95)
sigma_bayes = N_s_bayes / (luminosity_pb * signal_eff)

# Invert to get FM2 exclusion
x2_bayes = (sigma_bayes - b) / a
fM2_bayes = np.sqrt(x2_bayes) if x2_bayes > 0 else None

# === 2. Asimov 5σ Discovery (for comparison) ===
def discovery_significance_equation(s, B, z=5.0):
    return s / np.sqrt(s + B) - z

N_s_disc = fsolve(discovery_significance_equation, x0=5*np.sqrt(N_b), args=(N_b))[0]
sigma_disc = N_s_disc / (luminosity_pb * signal_eff)
x2_disc = (sigma_disc - b) / a
fM2_disc = np.sqrt(x2_disc) if x2_disc > 0 else None

# === Plotting ===
fM2_range = np.linspace(0, 1000, 500)
sigma_fM2 = a * fM2_range**2 + b

plt.figure(figsize=(8, 5))
plt.plot(fM2_range, sigma_fM2, label=r"$\sigma(FM2) = a x^2 + b$")
plt.axhline(sigma_bayes, color='blue', linestyle='--', label=r"Bayesian 95% CL Excl. $\sigma$")
plt.axhline(sigma_disc, color='orange', linestyle='--', label=r"5$\sigma$ Discovery $\sigma$")
plt.xlabel(r"$|FM2|$ [$\mathrm{TeV}^{-4}$]")
plt.ylabel(r"Cross-section [pb]")
plt.title("Cross-Section vs |FM2| with Bayesian 95% and 5σ Discovery")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig("fm2_bayesian_vs_discovery.pdf")
plt.show()

# === Summary Output ===
print(f"Bayesian 95% CL σ_s upper limit: {sigma_bayes:.6f} pb")
if fM2_bayes:
    print(f"✅ Bayesian 95% CL exclusion: |FM2| < {fM2_bayes:.1f} TeV⁻⁴")
else:
    print("❌ No Bayesian exclusion — cross-section too small")

if fM2_disc:
    print(f"✅ 5σ Discovery threshold: |FM2| > {fM2_disc:.1f} TeV⁻⁴")
else:
    print("❌ No discovery threshold — signal too small")
