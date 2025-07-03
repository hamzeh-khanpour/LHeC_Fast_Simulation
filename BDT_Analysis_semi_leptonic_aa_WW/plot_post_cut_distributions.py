import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.gridspec as gridspec

import os


# Load post-cut data
signal_df = pd.read_csv("signal_after_cut_FM2.csv")
background_df = pd.read_csv("background_after_cut_FM2.csv")


# Directory to save plots
output_dir = "post_cut_distributions"
os.makedirs(output_dir, exist_ok=True)


# Features to compare
features = [
    "lepton_pt",
    "lepton_eta",
    "leading_jet_pt",
    "leading_jet_eta",
    "subleading_jet_eta",
    "missing_et",
    "delta_r",
    "jet_centrality",
    "delta_eta_jj",
    "m_w_leptonic",
    "m_w_hadronic",
    "pt_w_leptonic",
    "pt_w_hadronic",
    "delta_phi_lep_met",
    "mt_w_leptonic",
    "ht_total",
    "delta_phi_jj",
    "delta_phi_wl_wh",
    "delta_eta_wl_wh",
    "m_jj",
    "m_lvjj"
]


# Create PDF summary report with S/B panel
with PdfPages("post_cut_summary_FM2.pdf") as pdf:
    for feature in features:
        fig = plt.figure(figsize=(8, 8))
        gs = gridspec.GridSpec(2, 1, height_ratios=[3, 1], hspace=0.05)
        ax0 = plt.subplot(gs[0])
        ax1 = plt.subplot(gs[1], sharex=ax0)

        # Define binning based on combined range
        min_val = min(signal_df[feature].min(), background_df[feature].min())
        max_val = max(signal_df[feature].max(), background_df[feature].max())
        bins = 50
        bin_edges = np.linspace(min_val, max_val, bins + 1)
        bin_centers = 0.5 * (bin_edges[1:] + bin_edges[:-1])

        # Histogram values
        sig_vals, _ = np.histogram(signal_df[feature], bins=bin_edges, weights=signal_df["weight"])
        bkg_vals, _ = np.histogram(background_df[feature], bins=bin_edges, weights=background_df["weight"])

        # Top: signal & background
        ax0.hist(bin_edges[:-1], bins=bin_edges, weights=sig_vals,
                 label="Signal", histtype='stepfilled', edgecolor='crimson', alpha=0.6)
        ax0.hist(bin_edges[:-1], bins=bin_edges, weights=bkg_vals,
                 label="Background", histtype='stepfilled', edgecolor='navy', alpha=0.6)

        ax0.set_ylabel("Normalized Entries")
        ax0.set_title(f"{feature.replace('_', ' ').title()} After BDT Cut")
        ax0.legend()
        ax0.grid(True, linestyle="--", alpha=0.5)

        # Bottom: S / √B ratio
        ratio = np.divide(sig_vals, np.sqrt(bkg_vals + 1e-6))  # S / √B
        ax1.plot(bin_centers, ratio, drawstyle="steps-mid", color="black")
        ax1.set_ylabel("S / √B")
        ax1.set_xlabel(feature.replace("_", " ").title())
        ax1.grid(True, linestyle="--", alpha=0.5)

        plt.tight_layout()
        pdf.savefig()
        plt.savefig(os.path.join(output_dir, f"post_cut_{feature}_with_ratio_FM2.pdf"))
        plt.close()

    # Optional: Add metadata to PDF
    d = pdf.infodict()
    d['Title'] = 'Post-Cut Distribution Summary with Ratio Panels'
    d['Author'] = 'Your Analysis Pipeline'
    d['Subject'] = 'Signal vs Background and S / √B Distributions after BDT Cut'
    d['Keywords'] = 'BDT Cut Physics ML XGBoost LHeC'
    d['CreationDate'] = pd.Timestamp.today()

print("📄 Summary PDF with S / √B panels saved as post_cut_summary.pdf ✅")
