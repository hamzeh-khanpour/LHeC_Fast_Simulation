import ROOT
import numpy as np
import matplotlib.pyplot as plt

# Open ROOT file
file = ROOT.TFile.Open("output_histograms.root")

# Load signal histogram
signal_hist = file.Get("signal_FM2_Lambda4/hist_jet_centrality_FM2_Lambda4")

# Define backgrounds (match your ROOT keys)
backgrounds = {
    "aa_ww": "γγ → WW",
    "aa_ttbar": "γγ → tt̄ (×100)",
    "aa_tautau": "γγ → ττ",
    "aa_mumu": "γγ → μμ",
    "inclusive_ttbar": "Inclusive tt̄",
    "single_top": "Single Top",
    "w_production": "W Production",
    "z_production": "Z Production",
    "wwj": "WWj",
    "zzj": "ZZj (×100)",
    "wzj": "WZj"
}

# Plot setup
plt.figure(figsize=(10, 7))
plt.title("Jet Centrality: Signal vs Backgrounds")
plt.xlabel(r"Jet Centrality $C_{\mathrm{jets}} = \frac{|\eta_1 + \eta_2|}{2}$")
plt.ylabel("Normalized Entries")
plt.yscale("log")
plt.grid(True, linestyle="--", alpha=0.5)

# Plot signal
n_bins = signal_hist.GetNbinsX()
bin_edges = np.array([signal_hist.GetBinLowEdge(i+1) for i in range(n_bins+1)])
signal_values = np.array([signal_hist.GetBinContent(i+1) for i in range(n_bins)])
if signal_values.sum() > 0:
    signal_values /= signal_values.sum()
plt.step(bin_edges[:-1], signal_values, where='mid', label="Signal (FM2)", color='red', linewidth=2)

# Plot each background
colors = plt.cm.get_cmap("tab20", len(backgrounds))
for i, (key, label) in enumerate(backgrounds.items()):
    hist_name = f"hist_jet_centrality_{key}"
    hist = file.Get(f"background_{key}/{hist_name}")
    if hist:
        values = np.array([hist.GetBinContent(j+1) for j in range(n_bins)])
        if values.sum() > 0:
            values /= values.sum()
            plt.step(bin_edges[:-1], values, where='mid',
                     label=label, linewidth=1.5, color=colors(i))

plt.legend(fontsize=8, loc="upper right")
plt.tight_layout()
plt.savefig("jet_centrality_comparison_all.pdf", dpi=300)
plt.show()
