import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import xgboost as xgb
import ROOT

#-------------------------
# Settings
#-------------------------
n_generated = 1_000_000
luminosity_fb = 100.0
luminosity_pb = luminosity_fb * 1000.0  # [pb^-1]

signal_cross_section_pb = 0.014288200000000001    # [pb] FM2
sigma_SM = 0.0099465                              # [pb] (Standard Model-only signal prediction)
#  0.0099465   for  aa_ww_semi_leptonic_SM_NP_1_FMi_0
#  0.0150743   for  aa_ww_semi_leptonic_SM

# Cross sections (pb) for backgrounds
background_cross_sections_pb = {
    "aa_ww": 0.0099465,
#  0.0099465   for  aa_ww_semi_leptonic_SM_NP_1_FMi_0
#  0.0150743   for  aa_ww_semi_leptonic_SM
    "aa_ttbar": 4.824851e-03 / 100.0,
    "aa_tautau": 2.51510000,
    "aa_mumu": 2.57270000,
    "inclusive_ttbar": 0.0065764,
    "single_top": 1.36209,
    "w_production": 1.910288,
    "z_production": 0.24064758729900002,
    "wwj": 0.016080595320336195,
    "zzj": 6.694889944457796e-03 / 100.0,
    "wzj": 0.0023785292894910495
}

#-------------------------
# Load ML dataset
#-------------------------
df = pd.read_csv("ml_input_from_histograms.csv")
has_process = "process" in df.columns

X = df.drop(columns=["label", "weight"] + (["process"] if has_process else []))
y = df["label"]
weights = df["weight"]

#-------------------------
# Train BDT
#-------------------------
model = xgb.XGBClassifier(
    n_estimators=200,
    max_depth=4,
    learning_rate=0.05,
    use_label_encoder=False,
    eval_metric="logloss"
)
model.fit(X, y, sample_weight=weights)
df["bdt_score"] = model.predict_proba(X)[:, 1]

#-------------------------
# Optimal threshold
#-------------------------
thresholds = np.linspace(0.0, 1.0, 100)
best_threshold = max(
    thresholds,
    key=lambda t: df[df["bdt_score"] > t].query("label == 1")["weight"].sum()
    / np.sqrt(df[df["bdt_score"] > t]["weight"].sum() + 1e-9)
)
df_cut = df[df["bdt_score"] > best_threshold]

#-------------------------
# Signal stats
#-------------------------
s_total = df[df["label"] == 1]["weight"].sum()
s_after = df_cut[df_cut["label"] == 1]["weight"].sum()
efficiency_bdt = s_after / s_total if s_total > 0 else 0.0

# Signal preselection efficiency
f = ROOT.TFile.Open("output_histograms.root")
hist_preselected_sig = f.Get("signal_FM2_Lambda4/hist_jet_centrality_FM2_Lambda4")
n_preselected_sig = hist_preselected_sig.Integral() if hist_preselected_sig else 0
eff_preselection_sig = n_preselected_sig / n_generated
eff_total_sig = eff_preselection_sig * efficiency_bdt

#-------------------------
# Background stats
#-------------------------
df_bkg = df[df["label"] == 0].copy()
df_bkg_cut = df_cut[df_cut["label"] == 0].copy()
if has_process:
    df_bkg["process"] = df["process"]
    df_bkg_cut["process"] = df_cut["process"]
else:
    df_bkg["process"] = "Unknown"
    df_bkg_cut["process"] = "Unknown"

b_total = df_bkg_cut["weight"].sum()
bkg_cutflow = df_bkg_cut.groupby("process")["weight"].sum().sort_values(ascending=False)

# MC background efficiency (total)
eff_bkg_total = df_bkg_cut.shape[0] / df_bkg.shape[0]

#-------------------------
# Background efficiencies and yields by process
#-------------------------
print("\n🧾 Background efficiencies and projected yields (100 fb⁻¹):")
total_yield_manual = 0
for proc in df_bkg["process"].unique():
    suffix = "_production" if proc in {"wwj", "zzj", "wzj"} else ""
    hist_name = f"hist_jet_centrality_{proc}{suffix}"
    proc_hist = f.Get(f"background_{proc}/{hist_name}")
    n_pre_bkg = proc_hist.Integral() if proc_hist else 0
    n_cut_bkg = df_bkg_cut[df_bkg_cut["process"] == proc].shape[0]
    eff_pre_bkg = n_pre_bkg / n_generated
    eff_bdt_bkg = n_cut_bkg / n_pre_bkg if n_pre_bkg > 0 else 0.0
    eff_total_bkg = eff_pre_bkg * eff_bdt_bkg
    xs = background_cross_sections_pb.get(proc, 0.0)
    expected_yield = eff_total_bkg * xs * luminosity_fb * 1000.0
    total_yield_manual += expected_yield
    print(f"  {proc:<16}  pre: {eff_pre_bkg:.4f},  BDT: {eff_bdt_bkg:.4f},  total: {eff_total_bkg:.4f},  yield: {expected_yield:.2f}")

f.Close()

#-------------------------
# Report
#-------------------------
print(f"\n💎 Optimal BDT cut: {best_threshold:.3f}")
print(f"🧬 BDT efficiency (from preselected): {efficiency_bdt:.4f}")
print(f"🧬 Preselection efficiency (from ROOT): {eff_preselection_sig:.4f} ({int(n_preselected_sig)} / {n_generated})")
print(f"🧬 Total signal efficiency (MC from 1M): {eff_total_sig:.4f}")
print(f"🎭 MC background efficiency (from all samples): {eff_bkg_total:.4f}")
print(f"📈 Significance Z = S / √(S+B): {s_after / np.sqrt(s_after + b_total):.2f}")
print(f"🎯 Total background after cut: {b_total:.2f} (weighted sum)")
print(f"📊 Total background yield from per-process calculation: {total_yield_manual:.2f}")
print("\n📊 Background yields by process (from DataFrame sum):")
print(bkg_cutflow.to_string())

#-------------------------
# Save filtered events
#-------------------------
df_cut[df_cut["label"] == 1].to_csv("signal_after_cut.csv", index=False)
df_cut[df_cut["label"] == 0].to_csv("background_after_cut.csv", index=False)

#-------------------------
# Plot BDT score
#-------------------------
plt.figure(figsize=(8, 6))
plt.hist(df[df["label"] == 1]["bdt_score"], bins=50, alpha=0.5, label="Signal", density=True)
plt.hist(df[df["label"] == 0]["bdt_score"], bins=50, alpha=0.5, label="Background", density=True)
plt.axvline(best_threshold, color='red', linestyle='--', label=f"Cut = {best_threshold:.3f}")
plt.xlabel("BDT Score")
plt.ylabel("Normalized Entries")
plt.title("BDT Score with Optimal Cut")
plt.legend()
plt.grid(True, linestyle="--", alpha=0.5)
plt.tight_layout()
plt.savefig("bdt_cut_applied.pdf")
plt.show()


