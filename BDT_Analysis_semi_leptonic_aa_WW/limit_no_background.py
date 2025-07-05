# ================================================================================
#   Hamzeh Khanpour — June 2025
#   No-Background 95% CL Limit on FM2 (cut-and-count + profile likelihood)
# ================================================================================

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import poisson
from scipy.optimize import root_scalar

# -----------------------------
# Inputs (Update if needed)
# -----------------------------
luminosity_fb = 1000.0  # Integrated luminosity in fb^-1
eff_sig = 0.1405        # Total signal efficiency = preselection × ML
sigma_sig_fb = 0.0173508 * 1000.0  # SM signal cross-section (FM2 = 0) [fb]
s95_bayesian = 3.0       # Bayesian Poisson upper limit at 95% CL, assuming n_obs = 0

# -----------------------------
# Cross section as function of FM2 [TeV^-4]
# -----------------------------
def sigma_fm2_fb(fm2):
    a = -1.318084e-02
    b = 7.832913e-04
    c = 1.491861e+01
    return a * fm2 + b * fm2**2 + c  # [fb]

# -----------------------------
# Compute signal yield for a given FM2
# -----------------------------
def get_signal_events(fm2):
    return sigma_fm2_fb(fm2) * eff_sig * luminosity_fb

# -----------------------------
# Profile likelihood definition
# -----------------------------
def q_mu(fm2):
    s = get_signal_events(fm2)
    b = 0.0  # no-background assumption
    n_obs = 0  # observed events
    lam = s + b
    if lam <= 0:
        return np.inf
    return -2.0 * (poisson.logpmf(n_obs, lam) - poisson.logpmf(n_obs, b + 1e-6))

# -----------------------------
# Attempt profile likelihood limit
# -----------------------------
def scan_profile_likelihood():
    print("\n🔎 Attempting profile likelihood limit...")
    def q_target(fm2): return q_mu(fm2) - 3.84
    try:
        upper = root_scalar(q_target, bracket=(0.01, 50), method="brentq").root
        print(f"✅ Profile Likelihood 95% CL Limit on FM2: fM2 < {upper:.4f} [TeV⁻⁴]")
        return upper
    except Exception as e:
        print("❌ Profile likelihood failed:", e)
        return None

# -----------------------------
# Bayesian inversion
# -----------------------------
def invert_sigma(s_target):
    def func(fm2):
        return sigma_fm2_fb(fm2) * eff_sig * luminosity_fb - s_target
    try:
        result = root_scalar(func, bracket=(0.001, 500), method="brentq")
        return result.root if result.converged else None
    except Exception as e:
        print("⚠️ Bayesian inversion failed:", e)
        return None

# -----------------------------
# Run
# -----------------------------
print("\n===================================")
print("      No-Background Limit Setting")
print("===================================")
print(f"SM Background Yield (aa→WW): {0.0149219 * 1000 * 0.1405 * 1000:.2f}")
print(f"Signal Efficiency × Luminosity: {eff_sig:.4f} × {luminosity_fb} fb⁻¹")

# Try profile likelihood
fm2_limit_profile = scan_profile_likelihood()

# Try Bayesian cut-and-count if needed
fm2_limit_bayes = invert_sigma(s95_bayesian)
if fm2_limit_bayes:
    print(f"📏 Bayesian 95% CL Upper Limit on FM2: fM2 < {fm2_limit_bayes:.4f} [TeV⁻⁴]")

# -----------------------------
# Plot likelihood scan
# -----------------------------
fm2_vals = np.linspace(0.01, 50, 200)
q_vals = [q_mu(f) for f in fm2_vals]

plt.figure(figsize=(8, 6))
plt.plot(fm2_vals, q_vals, label=r"$q(f_{M2}/\Lambda^4)$")
plt.axhline(3.84, color='red', linestyle='--', label="95% CL threshold")
if fm2_limit_profile:
    plt.axvline(fm2_limit_profile, color='gray', linestyle=':', label=f'Limit: {fm2_limit_profile:.2f}')
plt.xlabel(r"$f_{M2}/\Lambda^4$ [$\mathrm{TeV}^{-4}$]")
plt.ylabel(r"$q(f_{M2}/\Lambda^4)$")
plt.title("Profile Likelihood Scan (No Background)")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig("limit_no_background_likelihood_scan.pdf")
plt.show()

