import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import poisson
from scipy.optimize import root_scalar

# ================================
# 🎯 Fixed Parameters
# ================================
luminosity_fb = 100.0
signal_eff = 0.2054
background_eff = 0.1692
sigma_background_fb = 9.9465
delta_sys_frac = 0.10
n_bins = 20

# ===================================
# 📥 Load Preprocessed Data (Mocked)
# ===================================
np.random.seed(42)
n_events = 5000
data = pd.DataFrame({
    "bdt_score": np.random.uniform(0, 1, n_events),
    "label": np.random.choice([0, 1], size=n_events, p=[0.8, 0.2]),
    "weight": np.random.exponential(scale=1.0, size=n_events)
})
signal_template = data[data["label"] == 1]
background = data[data["label"] == 0]

# Histogram binning
bin_edges = np.linspace(0, 1, n_bins + 1)
b_hist, _ = np.histogram(background["bdt_score"], bins=bin_edges, weights=background["weight"])
b_hist_scaled = b_hist * (sigma_background_fb * background_eff * luminosity_fb / np.sum(background["weight"]))
n_obs = np.round(b_hist_scaled)  # Asimov dataset (background-only hypothesis)

# ===================================
# 📐 σ(fM2) function in femtobarns
# ===================================
def sigma_fm2_fb(fm2):
    return -8.777e-3 * fm2 + 0.5215 * fm2**2 + 9.941

# Use normalized shape of signal BDT distribution as template
s_hist_template, _ = np.histogram(signal_template["bdt_score"], bins=bin_edges, weights=signal_template["weight"])
s_hist_template = s_hist_template / np.sum(signal_template["weight"])  # normalized shape

# ===================================
# 🔬 Likelihood Definitions
# ===================================
def get_scaled_signal_hist(fm2):
    sigma = sigma_fm2_fb(fm2)
    scale = sigma * signal_eff * luminosity_fb
    return s_hist_template * scale

def likelihood(fm2):
    s = get_scaled_signal_hist(fm2)
    L = 0.0
    for i in range(n_bins):
        lam = s[i] + b_hist_scaled[i]
        L += poisson.logpmf(n_obs[i], lam)
    return np.exp(L)

def q_mu(fm2):
    L_mu = likelihood(fm2)
    L_SM = likelihood(0.0)
    if L_mu == 0 or L_SM == 0:
        return np.inf
    return -2.0 * np.log(L_mu / L_SM)

# =============================
# 🔎 Find 95% CL Interval
# =============================
def q_target(fm2):
    return q_mu(fm2) - 3.84  # 95% CL threshold

try:
    fm2_upper = root_scalar(q_target, bracket=(0.0001, 2.0), method='brentq').root
    fm2_lower = root_scalar(q_target, bracket=(-2.0, -0.0001), method='brentq').root
except:
    fm2_upper = None
    fm2_lower = None

# =============================
# 📈 Plot q(FM2)
# =============================
fm2_vals = np.linspace(-1.0, 1.0, 200)
q_vals = [q_mu(fm2) for fm2 in fm2_vals]

plt.figure(figsize=(8, 5))
plt.plot(fm2_vals, q_vals, label=r'$q(f_{M2})$')
plt.axhline(3.84, color='red', linestyle='--', label='95% CL Threshold')
if fm2_lower and fm2_upper:
    plt.axvline(fm2_lower, color='gray', linestyle=':', label=f'Lower limit: {fm2_lower:.4f}')
    plt.axvline(fm2_upper, color='gray', linestyle=':', label=f'Upper limit: {fm2_upper:.4f}')
plt.xlabel(r'$f_{M2}$ [$\mathrm{TeV}^{-4}$]')
plt.ylabel(r'$q(f_{M2})$')
plt.title('Profile Likelihood Scan over $f_{M2}$')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()

# =============================
# 🖨️ Summary
# =============================
print("\n======== Shape-Based FM2 Limit Scan ========")
if fm2_lower and fm2_upper:
    print(f"📉 95% CL Excluded FM2 Range: {fm2_lower:.4f} to {fm2_upper:.4f} [TeV^-4]")
else:
    print("❌ No valid FM2 limits found in the scanned range.")
