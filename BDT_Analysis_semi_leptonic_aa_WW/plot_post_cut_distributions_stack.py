import pandas as pd
import numpy as np  # ✅ Missing import
import matplotlib.pyplot as plt

# Load post-cut data
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

    # Define binning based on the combined range
    min_val = min(signal_df[feature].min(), background_df[feature].min())
    max_val = max(signal_df[feature].max(), background_df[feature].max())
    bins = 50
    bin_edges = np.linspace(min_val, max_val, bins + 1)

    # Histogram values
    sig_vals, _ = np.histogram(signal_df[feature], bins=bin_edges, weights=signal_df["weight"])
    bkg_vals, _ = np.histogram(background_df[feature], bins=bin_edges, weights=background_df["weight"])

    # Stack bars
    plt.hist([bin_edges[:-1], bin_edges[:-1]], bins=bin_edges,
             weights=[bkg_vals, sig_vals],
             label=["Background", "Signal"],
             stacked=True, color=["steelblue", "crimson"],
             histtype='stepfilled', alpha=0.8)

    plt.xlabel(feature.replace("_", " ").title())
    plt.ylabel("Weighted Event Yield")
    plt.title(f"Stacked Post-Cut Distribution: {feature.replace('_', ' ').title()}")
    plt.legend()
    plt.grid(True, linestyle="--", alpha=0.5)
    plt.tight_layout()
    plt.savefig(f"post_cut_{feature}_stacked.pdf")
    plt.close()

print("✅ All stacked post-cut distributions saved as PDF.")
