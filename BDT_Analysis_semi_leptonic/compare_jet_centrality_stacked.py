import ROOT
import numpy as np
import matplotlib.pyplot as plt

# Open ROOT file
file = ROOT.TFile.Open("output_histograms.root")

# Signal
signal_hist = file.Get("signal_FM2_Lambda4/hist_jet_centrality_FM2_Lambda4")

# Backgrounds (match with ROOT structure)
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

# Handle suffixes for special names
special_suffix = {
    "wwj": "_production",
    "zzj": "_production",
    "wzj": "_production"
}

# Bin info
n_bins = signal_hist.GetNbinsX()
bin_edges = np.array([signal_hist.GetBinLowEdge(i+1) for i in range(n_bins+1)])
bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
bin_width = bin_edges[1] - bin_edges[0]

# Collect background stacks
bkg_labels = []
bkg_stacks = []
colors = plt.cm.get_cmap("tab20", len(backgrounds))

for i, (key, label) in enumerate(backgrounds.items()):
    suffix = special_suffix.get(key, "")
    hist_path = f"background_{key}/hist_jet_centrality_{key}{suffix}"
    hist = file.Get(hist_path)
    if hist:
        values = np.array([hist.GetBinContent(j+1) for j in range(n_bins)])
        if values.sum() > 0:
            bkg_labels.append(label)
            bkg_stacks.append(values)
        else:
            print(f"⚠️ Histogram for {key} is empty.")
    else:
        print(f"❌ Histogram not found: {hist_path}")

# Convert background stack
bkg_stacks = np.array(bkg_stacks)
stack_bottom = np.zeros(n_bins)

# Plotting
plt.figure(figsize=(10, 7))
plt.title("Jet Centrality: Stacked Backgrounds + Signal Overlay")
plt.xlabel(r"Jet Centrality $C_{\mathrm{jets}} = \frac{|\eta_1 + \eta_2|}{2}$")
plt.ylabel("Events")
plt.yscale("log")
plt.grid(True, linestyle="--", alpha=0.5)

# Plot backgrounds
for i, label in enumerate(bkg_labels):
    plt.bar(bin_centers, bkg_stacks[i], width=bin_width, bottom=stack_bottom,
            label=label, color=colors(i), alpha=0.8, edgecolor='black', linewidth=0.3)
    stack_bottom += bkg_stacks[i]

# Overlay signal, scaled to background yield
signal_values = np.array([signal_hist.GetBinContent(i+1) for i in range(n_bins)])
if signal_values.sum() > 0:
    scale = stack_bottom.sum() / signal_values.sum()
    scaled_signal = signal_values * scale
    plt.step(bin_edges[:-1], scaled_signal, where="mid", color="black",
             linewidth=2, label="Signal (FM2)")

# Final touches
plt.legend(fontsize=8, loc="upper right")
plt.tight_layout()
plt.savefig("jet_centrality_stacked_all.pdf", dpi=300)
plt.show()
