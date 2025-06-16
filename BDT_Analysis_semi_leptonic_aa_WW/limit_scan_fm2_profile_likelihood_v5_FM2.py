

## ================================================================================
##        Hamzeh Khanpour  --- June 2025
## ================================================================================



import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import poisson
from scipy.optimize import root_scalar
import ROOT


# ================================
# 🎯 Fixed Parameters
# ================================
luminosity_fb = 100.0

sigma_background_fb = 9.9465  # SM-only gamma-gamma -> WW [fb]

delta_sys_frac = 0.10
n_bins = 20
n_generated = 1_000_000  # Total generated events



# Background preselection efficiencies
background_cross_sections_fb = {
    "aa_ww": 0.0099465 * 1000.0,
    "aa_ttbar": 4.824851e-03 / 100.0 * 1000.0,
    "aa_tautau": 2.51510000 * 1000.0,
    "aa_mumu": 2.57270000 * 1000.0,
    "inclusive_ttbar": 0.0065764 * 1000.0,
    "single_top": 1.36209 * 1000.0,
    "w_production": 1.910288 * 1000.0,
    "z_production": 0.24064758729900002 * 1000.0,
    "wwj": 0.016080595320336195 * 1000.0,
    "zzj": 6.694889944457796e-03 / 100.0 * 1000.0,
    "wzj": 0.0023785292894910495 * 1000.0
}



# ===================================
# 📐 Cross section function σ(fM2)
# ===================================
def sigma_fm2_fb(fm2):
    a =  -8.615e-3
    b =  5.221e-4
    c =  9.9465
    return a * fm2 + b * fm2**2 + c



# ================================
# 📥 Load Data
# ================================
df = pd.read_csv("ml_with_bdt_scores_FM2.csv")
has_process = "process" in df.columns



# =================================
# 🔍 Optimize BDT Cut
# =================================
bdt_cuts = np.linspace(0.0, 1.0, 100)
ams_scores = []

for cut in bdt_cuts:
    s = df[(df["label"] == 1) & (df["bdt_score"] > cut)]["weight"].sum()
    b = df[(df["label"] == 0) & (df["bdt_score"] > cut)]["weight"].sum()
    ams = s / np.sqrt(b + 1e-6) if b > 0 else 0.0
    ams_scores.append(ams)

idx_best = np.argmax(ams_scores)
bdt_cut = bdt_cuts[idx_best]
print(f"✅ Optimal BDT Cut: {bdt_cut:.3f}")



# =================================
# 🔪 Apply BDT Cut
# =================================
df_cut = df[df["bdt_score"] > bdt_cut]
signal_cut = df_cut[df_cut["label"] == 1]
background_cut = df_cut[df_cut["label"] == 0]
signal = df[df["label"] == 1]
background = df[df["label"] == 0]



# BDT efficiencies
signal_eff_ml = len(signal_cut) / len(signal)
background_eff_ml = len(background_cut) / len(background)

print(f"✅ ML Efficiencies → signal: {signal_eff_ml:.4f}, background: {background_eff_ml:.4f}")



# =================================
# ⚙️ Preselection Efficiencies
# =================================
root_file = ROOT.TFile.Open("output_histograms_FM2.root")

# Signal preselection efficiency
hist_sig_pre = root_file.Get("signal_FM2_Lambda4/hist_jet_centrality_FM2_Lambda4")
n_sig_pre = hist_sig_pre.Integral() if hist_sig_pre else 0
eff_preselection_sig = n_sig_pre / n_generated
signal_eff_total = eff_preselection_sig * signal_eff_ml
print(f"✅ Total Signal Efficiency (Preselection × ML): {signal_eff_total:.4f}")



preselection_counts = {}
for bkg in background_cross_sections_fb:
    path = f"background_{bkg}/hist_jet_centrality_{bkg}"
    hist = root_file.Get(path)
    n_pre = hist.Integral() if hist else 0
    preselection_counts[bkg] = n_pre

eff_preselection_bkg = {
    bkg: preselection_counts[bkg] / n_generated
    for bkg in preselection_counts
}



total_preselection_bkg_eff = 0.0
for bkg in background_cut["process"].unique():
    frac = (background_cut["process"] == bkg).sum() / len(background_cut)
    total_preselection_bkg_eff += frac * eff_preselection_bkg.get(bkg, 1.0)

background_eff_total = background_eff_ml * total_preselection_bkg_eff
print(f"✅ Total Background Efficiency (Preselection × ML): {background_eff_total:.4f}")





# ===================================
# 📊 Build Histograms
# ===================================
bin_edges = np.linspace(0, 1, n_bins + 1)
s_hist_template, _ = np.histogram(signal_cut["bdt_score"], bins=bin_edges, weights=signal_cut["weight"])
s_hist_template = s_hist_template / np.sum(s_hist_template)

b_hist, _ = np.histogram(background_cut["bdt_score"], bins=bin_edges, weights=background_cut["weight"])
b_total_weight = np.sum(background_cut["weight"])
b_hist_scaled = b_hist * (sigma_background_fb * background_eff_total * luminosity_fb / b_total_weight)
n_obs = np.round(b_hist_scaled)  # Asimov background-only



# ===================================
# 🔬 Likelihood Definitions
# ===================================
def get_scaled_signal_hist(fm2):
    sigma = sigma_fm2_fb(fm2)
    scale = sigma * signal_eff_total * luminosity_fb
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



# ===================================
# 🔎 Find 95% CL Interval
# ===================================
def q_target(fm2):
    return q_mu(fm2) - 3.84

try:
    fm2_upper = root_scalar(q_target, bracket=(0.001, 50.0), method='brentq').root
    fm2_lower = root_scalar(q_target, bracket=(-50.0, -0.001), method='brentq').root
except:
    fm2_upper = None
    fm2_lower = None



# ===================================
# 📈 Plot q(FM2)
# ===================================
fm2_vals = np.linspace(-50.0, 50.0, 200)
q_vals = [q_mu(fm2) for fm2 in fm2_vals]

plt.figure(figsize=(8, 6))

plt.plot(fm2_vals, q_vals, label=r'$q(f_{M1}/\Lambda^4)$', linewidth=2)
plt.axhline(3.84, color='red', linestyle='--', label='95% CL Threshold', linewidth=2)
if fm2_lower and fm2_upper:
    plt.axvline(fm2_lower, color='gray', linestyle=':', label=f'Lower limit: {fm2_lower:.4f}', linewidth=2)
    plt.axvline(fm2_upper, color='gray', linestyle=':', label=f'Upper limit: {fm2_upper:.4f}', linewidth=2)

plt.xlabel(r'$f_{M2}/\Lambda^4$ [$\mathrm{TeV}^{-4}$]')
plt.ylabel(r'$q(f_{M2}/\Lambda^4)$')
plt.title('Profile Likelihood Scan over $f_{M2}/\Lambda^4$')
plt.legend()
plt.grid(True)
plt.tight_layout()

plt.savefig("limit_scan_fm2_profile_likelihood_v5_FM2.pdf")
print("✅ Saved as 'limit_scan_fm2_profile_likelihood_v5_FM2.pdf'")
plt.show()




# ===================================
# 🖨️ Final Report
# ===================================
print("\n======== Final FM1 Limit Scan Report ========")
print(f"Signal Efficiency (total):    {signal_eff_total:.4f}")
print(f"Background Efficiency (total): {background_eff_total:.4f}")
if fm2_lower and fm2_upper:
    print(f"📉 95% CL Excluded FM2 Range: {fm2_lower:.4f} to {fm2_upper:.4f} [TeV^-4]")
else:
    print("❌ No valid FM1 limits found in the scanned range.")


