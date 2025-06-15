import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import poisson
from scipy.optimize import root_scalar



# ================================
# 📥 Load Data
# ================================
df = pd.read_csv("ml_with_bdt_scores.csv")
signal = df[df["label"] == 1]
background = df[df["label"] == 0]



# ================================
# 🧠 Optimize BDT Cut
# ================================
cut_values = np.linspace(0, 1, 200)
best_significance = 0
best_cut = 0.0

for cut in cut_values:
    s_cut = signal[signal["bdt_score"] > cut]
    b_cut = background[background["bdt_score"] > cut]
    s = np.sum(s_cut["weight"])
    b = np.sum(b_cut["weight"])
    if b > 0:
        Z = s / np.sqrt(s + b)
        if Z > best_significance:
            best_significance = Z
            best_cut = cut

bdt_cut = best_cut
print(f"✅ Optimal BDT Cut = {bdt_cut:.3f} (Z = {best_significance:.2f})")



# ================================
# 🎯 Fixed Parameters
# ================================
luminosity_fb = 100.0
sigma_background_fb = 9.9465
delta_sys_frac = 0.10
n_bins = 20



# ================================
# 📐 σ(fM2) Function
# ================================
def sigma_fm2_fb(fm2):
    a = -8.777e-3
    b = 0.5215
    c = 9.9465
    return a * fm2 + b * fm2**2 + c



# ================================
# 🔪 Apply BDT Cut
# ================================
signal_cut = signal[signal["bdt_score"] > bdt_cut]
background_cut = background[background["bdt_score"] > bdt_cut]

signal_eff = len(signal_cut) / len(signal)
background_eff = len(background_cut) / len(background)
print(f"✅ Efficiencies → signal: {signal_eff:.4f}, background: {background_eff:.4f}")




# ================================
# 📊 Histogram Binning
# ================================
bin_edges = np.linspace(0, 1, n_bins + 1)
s_hist_template, _ = np.histogram(signal_cut["bdt_score"], bins=bin_edges, weights=signal_cut["weight"])
s_hist_template /= np.sum(s_hist_template)

b_hist, _ = np.histogram(background_cut["bdt_score"], bins=bin_edges, weights=background_cut["weight"])
b_total_weight = np.sum(background_cut["weight"])
b_hist_scaled = b_hist * (sigma_background_fb * background_eff * luminosity_fb / b_total_weight)
n_obs = np.round(b_hist_scaled)




# ================================
# 🔬 Likelihood
# ================================
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



# ================================
# 🔎 Find 95% CL Interval
# ================================
def q_target(fm2):
    return q_mu(fm2) - 3.84

try:
    fm2_upper = root_scalar(q_target, bracket=(0.01, 2.0), method='brentq').root
    fm2_lower = root_scalar(q_target, bracket=(-2.0, -0.01), method='brentq').root
except:
    fm2_upper = None
    fm2_lower = None



# ================================
# 📈 Plot q(FM2)
# ================================
fm2_vals = np.linspace(-1.0, 1.0, 200)
q_vals = [q_mu(fm2) for fm2 in fm2_vals]

plt.figure(figsize=(8, 5))
plt.plot(fm2_vals, q_vals, label=r'$q(f_{M2})$')
plt.axhline(3.84, color='red', linestyle='--', label='95% CL Threshold')
if fm2_lower and fm2_upper:
    plt.axvline(fm2_lower, color='gray', linestyle=':', label=f'{fm2_lower:.4f}')
    plt.axvline(fm2_upper, color='gray', linestyle=':', label=f'{fm2_upper:.4f}')
plt.xlabel(r'$f_{M2}$ [$\mathrm{TeV}^{-4}$]')
plt.ylabel(r'$q(f_{M2})$')
plt.title('Profile Likelihood Scan over $f_{M2}$')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig("limit_scan_fm2_profile_likelihood_v4_final.pdf")
plt.show()



# ================================
# 🖨️ Final Report
# ================================
print("\n======== Final FM2 Limit Scan Report ========")
print(f"Signal Efficiency:    {signal_eff:.4f}")
print(f"Background Efficiency: {background_eff:.4f}")
if fm2_lower and fm2_upper:
    print(f"📉 95% CL Excluded FM2 Range: {fm2_lower:.4f} to {fm2_upper:.4f} [TeV^-4]")
else:
    print("❌ No valid FM2 limits found in the scanned range.")



