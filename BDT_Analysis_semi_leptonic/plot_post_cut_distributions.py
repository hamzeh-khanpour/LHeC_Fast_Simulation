import pandas as pd
import matplotlib.pyplot as plt

# Load data
signal_df = pd.read_csv("signal_after_cut.csv")
background_df = pd.read_csv("background_after_cut.csv")

# Features to compare
features = [
    "lepton_pt", "lepton_eta",
    "leading_jet_pt", "leading_jet_eta", "subleading_jet_eta",
    "missing_et", "delta_r", "jet_centrality", "delta_eta_jj",
    "m_w_leptonic", "m_w_hadronic"
]

# Plot for each feature
for feature in features:
    plt.figure(figsize=(8, 6))
    plt.hist(signal_df[feature], bins=50, density=True, alpha=0.6,
             label="Signal", histtype='stepfilled', edgecolor='crimson')
    plt.hist(background_df[feature], bins=50, density=True, alpha=0.6,
             label="Background", histtype='stepfilled', edgecolor='navy')

    plt.xlabel(feature.replace("_", " ").title())
    plt.ylabel("Normalized Entries")
    plt.title(f"{feature.replace('_', ' ').title()} After BDT Cut")
    plt.legend()
    plt.grid(True, linestyle="--", alpha=0.5)
    plt.tight_layout()
    plt.savefig(f"post_cut_{feature}.pdf")
    plt.close()

print("✅ All post-cut distributions saved as PDF.")
