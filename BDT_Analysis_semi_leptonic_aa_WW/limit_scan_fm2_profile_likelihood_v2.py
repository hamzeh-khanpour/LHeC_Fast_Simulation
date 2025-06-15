import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import poisson
from scipy.optimize import root_scalar

# ================================
# 🎯 Fixed Parameters
# ================================
luminosity_fb = 100.0
signal_eff = 0.4828          # total signal efficiency from Delphes+BDT
background_eff = 0.0295      # total background efficiency from Delphes+BDT
sigma_background_fb = 9.9465
delta_sys_frac = 0.10
n_bins = 20



# ===================================
# 📥 Load Real Data
# ===================================
data = pd.read_csv("ml_with_bdt_scores.csv")
assert {'bdt_score', 'label', 'weight'}.issubset(data.columns), "CSV must have bdt_score, label, weight columns"


signal_template = data[data["label"] == 1]
background = data[data["label"] == 0]


# Histogram binning for BDT scores
bin_edges = np.linspace(0, 1, n_bins + 1)
b_hist, _ = np.histogram(background["bdt_score"], bins=bin_edges, weights=background["weight"])
b_hist_scaled = b_hist * (sigma_background_fb * background_eff * luminosity_fb / np.sum(background["weight"]))
n_obs = np.round(b_hist_scaled)  # Asimov dataset = background-only expectation



# ===================================
# 📐 σ(fM2) function in femtobarns
# ===================================
def sigma_fm2_fb(fm2):
    a = -8.777e-3        # fb / TeV^4
    b = 0.5215           # fb / TeV^8
    c = 9.9465           # fb
    return a * fm2 + b * fm2**2 + c



# ===================================
# 🧪 Signal Template (normalized shape from FM2 benchmark)
# ===================================
s_hist_template, _ = np.histogram(signal_template["bdt_score"], bins=bin_edges, weights=signal_template["weight"])
s_hist_template = s_hist_template / np.sum(s_hist_template)



# ===================================
# 🔬 Likelihood Definitions
# ===================================
def get_scaled_signal_hist(fm2):
    sigma = sigma_fm2_fb(fm2)
    return s_hist_template * sigma * signal_eff * luminosity_fb

def likelihood(fm2):
    s = get_scaled_signal_hist(fm2)
    logL = 0.0
    for i in range(n_bins):
        mu = s[i] + b_hist_scaled[i]
        logL += poisson.logpmf(n_obs[i], mu)
    return np.exp(logL)

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
    return q_mu(fm2) - 3.84  # 95% CL

try:
    fm2_upper = root_scalar(q_target, bracket=(0.01, 2.0), method='brentq').root
    fm2_lower = root_scalar(q_target, bracket=(-2.0, -0.01), method='brentq').root
except:
    fm2_upper, fm2_lower = None, None



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
plt.savefig("profile_likelihood_fm2_plot.pdf")
print("✅ Saved: profile_likelihood_fm2_plot.pdf")
plt.show()



# =============================
# 🖨️ Summary
# =============================
print("\n======== Final FM2 95% CL Exclusion Limits ========")
if fm2_lower and fm2_upper:
    print(f"📉 95% CL Excluded FM2 Range: {fm2_lower:.4f} to {fm2_upper:.4f} [TeV^-4]")
else:
    print("❌ No valid FM2 limits found in the scanned range.")


