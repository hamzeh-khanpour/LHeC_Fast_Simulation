
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import poisson, norm
from scipy.optimize import minimize_scalar, root_scalar

# ================================
# 🎯 Fixed Parameters
# ================================
luminosity_fb = 100.0
signal_eff = 0.2054
background_eff = 0.1692
sigma_background_fb = 9.9465
n_bins = 20
background_systematic_frac = 0.10

# ===================================
# 📥 Load Real BDT Scores
# ===================================
df = pd.read_csv("ml_input_from_histograms.csv")
signal = df[df["label"] == 1]
background = df[df["label"] == 0]

bin_edges = np.linspace(0, 1, n_bins + 1)
b_hist, _ = np.histogram(background["bdt_score"], bins=bin_edges, weights=background["weight"])
b_total = np.sum(background["weight"])
b_hist_scaled = b_hist * (sigma_background_fb * background_eff * luminosity_fb / b_total)
n_obs = np.round(b_hist_scaled)

# ===================================
# 📐 σ(fM2) function
# ===================================
def sigma_fm2_fb(fm2):
    return -8.777e-3 * fm2 + 0.5215 * fm2**2 + 9.941

s_hist_template, _ = np.histogram(signal["bdt_score"], bins=bin_edges, weights=signal["weight"])
s_hist_template = s_hist_template / np.sum(signal["weight"])

# ===================================
# 🔬 Nuisance-aware Likelihood
# ===================================
def get_scaled_signal_hist(fm2):
    sigma = sigma_fm2_fb(fm2)
    scale = sigma * signal_eff * luminosity_fb
    return s_hist_template * scale

def negative_log_likelihood(fm2, theta):
    s = get_scaled_signal_hist(fm2)
    b = b_hist_scaled * theta
    nll = 0.0
    for i in range(n_bins):
        mu_i = s[i] + b[i]
        if mu_i > 0:
            nll -= poisson.logpmf(n_obs[i], mu_i)
    # Gaussian constraint on theta (nuisance)
    nll -= norm.logpdf(theta, loc=1.0, scale=background_systematic_frac)
    return nll

# ===================================
# Profile Likelihood Ratio q(fM2)
# ===================================
def q_mu(fm2):
    min_nom = minimize_scalar(lambda t: negative_log_likelihood(fm2, t), bounds=(0.5, 1.5), method='bounded')
    min_den = minimize_scalar(lambda t: negative_log_likelihood(0.0, t), bounds=(0.5, 1.5), method='bounded')
    if not (min_nom.success and min_den.success):
        return np.inf
    return 2 * (min_nom.fun - min_den.fun)

# ===================================
# 95% CL scan
# ===================================
def q_target(fm2):
    return q_mu(fm2) - 3.84

try:
    fm2_upper = root_scalar(q_target, bracket=(0.01, 2.0), method='brentq').root
    fm2_lower = root_scalar(q_target, bracket=(-2.0, -0.01), method='brentq').root
except:
    fm2_upper = None
    fm2_lower = None

# ===================================
# Plot
# ===================================
fm2_vals = np.linspace(-0.5, 0.5, 200)
q_vals = [q_mu(f) for f in fm2_vals]

plt.figure(figsize=(8, 5))
plt.plot(fm2_vals, q_vals, label=r'$q(f_{M2})$ with systematics')
plt.axhline(3.84, color='red', linestyle='--', label='95% CL threshold')
if fm2_lower and fm2_upper:
    plt.axvline(fm2_lower, color='gray', linestyle=':', label=f'{fm2_lower:.4f}')
    plt.axvline(fm2_upper, color='gray', linestyle=':')
plt.xlabel(r'$f_{M2}$ [$\mathrm{TeV}^{-4}$]')
plt.ylabel(r'$q(f_{M2})$')
plt.title('Profile Likelihood with Background Systematics')
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()

# ===================================
# Result
# ===================================
print("\n======== Profile Likelihood Limits with Systematics ========")
if fm2_lower and fm2_upper:
    print(f"95% CL exclusion: f_M2 ∈ [{fm2_lower:.4f}, {fm2_upper:.4f}] TeV^-4")
else:
    print("❌ No valid FM2 exclusion range found.")
