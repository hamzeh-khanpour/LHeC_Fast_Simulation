import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

# Matplotlib configuration for publication-quality plots
import mplhep as hep
hep.style.use("CMS")

# ----------------------------
# Configuration
# ----------------------------
input_csv = "ml_input_from_histograms_FM2.csv"
output_dir = "sanity_check_plots"
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
# Load CSV data
# ----------------------------
print(f"📥 Loading: {input_csv}")
df = pd.read_csv(input_csv)

# ----------------------------
# Generate plots
# ----------------------------
for obs in observables:
    if obs not in df.columns:
        print(f"⚠️ Skipping missing column: {obs}")
        continue

    plt.figure(figsize=(11, 12))
    plt.subplots_adjust(left=0.15, right=0.95, bottom=0.12, top=0.95)

    sns.histplot(
        data=df, x=obs, hue="label", bins=50,
        stat="density", weights=df["weight"],  # ✅ Apply cross-section × lumi weights
        element="step", common_norm=False,
        palette={0: "tab:blue", 1: "tab:red"}
    )

    plt.title(f"Sanity Check (Weighted): {obs}")
    plt.xlabel(obs)
    plt.ylabel("Normalized Entries (weighted)")
    plt.legend(title="Label", labels=["Background", "Signal"])
    plt.tight_layout()

    save_path = os.path.join(output_dir, f"{obs}.png")
    plt.savefig(save_path)
    plt.close()
    print(f"✅ Saved: {save_path}")

print("\n✅ All weighted plots generated in:", output_dir)
