import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os


import mplhep as hep

# Use CMS plot style
hep.style.use("CMS")



# ----------------------------
# Configuration
# ----------------------------
input_csv = "ml_input_from_histograms_FM2.csv"
output_dir = "differential_distributions_all_backgrounds"
os.makedirs(output_dir, exist_ok=True)

# ----------------------------
# List of observables
# ----------------------------
observables = [
    "lepton_pt", "leading_jet_pt", "lepton_eta", "delta_r", "missing_et",
    "subleading_jet_eta", "leading_jet_eta", "jet_centrality", "delta_eta_jj",
    "m_w_hadronic", "m_w_leptonic", "pt_w_leptonic", "pt_w_hadronic",
    "delta_phi_lep_met", "mt_w_leptonic", "ht_total", "delta_phi_jj",
    "delta_phi_wl_wh", "delta_eta_wl_wh", "m_jj", "m_lvjj"
]

# ----------------------------
# Load Data
# ----------------------------
print(f"📥 Loading: {input_csv}")
df = pd.read_csv(input_csv)

# ----------------------------
# Separate signal and background
# ----------------------------
signal_df = df[df["label"] == 1]
background_df = df[df["label"] == 0]
background_processes = background_df["process"].unique()



# ----------------------------
# Generate one plot per observable
# ----------------------------
for obs in observables:
    if obs not in df.columns:
        print(f"⚠️ Skipping missing column: {obs}")
        continue



    plt.figure(figsize=(11, 12))
    plt.subplots_adjust(left=0.15, right=0.95, bottom=0.12, top=0.95)



    # Plot signal
    sns.histplot(
        data=signal_df, x=obs, weights=signal_df["weight"],
        bins=50, stat="density", element="step",
        label="FM2 Signal", color="blue", fill=False, linewidth=2
    )



    # Plot all backgrounds
    for bkg in background_processes:
        bkg_df = background_df[background_df["process"] == bkg]
        sns.histplot(
            data=bkg_df, x=obs, weights=bkg_df["weight"],
            bins=50, stat="density", element="step",
            label=bkg, fill=False, linewidth=1, alpha=0.7
        )



    plt.xlabel(obs)
    plt.ylabel(r"Normalized Differential Cross Section [1/GeV or 1/bin]")
    plt.title(f"Differential Distribution: {obs} (Signal vs Backgrounds)")
    plt.legend(fontsize=9, loc="upper right", frameon=True)
    plt.grid(True, linestyle="--", alpha=0.5)

    filename = f"{obs}_FM2_vs_backgrounds.png".replace("/", "_")
    plt.savefig(os.path.join(output_dir, filename))
    plt.close()
    print(f"✅ Saved: {filename}")



print("\n✅ All combined overlay plots generated in:", output_dir)
