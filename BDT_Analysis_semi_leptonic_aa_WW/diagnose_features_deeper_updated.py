import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# LaTeX-style labels for nicer plotting
feature_labels = {
    "lepton_pt": r"$p_T^\ell$",
    "leading_jet_pt": r"$p_T^{j_1}$",
    "subleading_jet_pt": r"$p_T^{j_2}$",
    "lepton_eta": r"$\eta^\ell$",
    "leading_jet_eta": r"$\eta^{j_1}$",
    "subleading_jet_eta": r"$\eta^{j_2}$",
    "missing_et": r"$E_T^{\mathrm{miss}}$",
    "delta_r": r"$\Delta R_{\ell j}$",
    "jet_centrality": r"${C^j}$",
    "delta_eta_jj": r"$\Delta\eta_{jj}$",
    "m_w_leptonic": r"$m_W^{\mathrm{lep}}$",
    "m_w_hadronic": r"$m_W^{\mathrm{had}}$",
    "pt_w_leptonic": r"$p_T^{W^{\mathrm{lep}}}$",
    "pt_w_hadronic": r"$p_T^{W^{\mathrm{had}}}$",
    "delta_phi_lep_met": r"$\Delta\phi_{\ell,\mathrm{MET}}$",
    "mt_w_leptonic": r"$m_T^{W^{\mathrm{lep}}}$",
    "ht_total": r"$H_T$",
    "delta_phi_jj": r"$\Delta\phi_{jj}$",
    "delta_phi_wl_wh": r"$\Delta\phi_{W^{\mathrm{lep}}, W^{\mathrm{had}}}$",
    "delta_eta_wl_wh": r"$\Delta\eta_{W^{\mathrm{lep}}, W^{\mathrm{had}}}$",
    "m_jj": r"$m_{jj}$",
    "m_lvvjj": r"$m_{{\ell \nu jj}}$"
}

# Load data
df = pd.read_csv("ml_input_from_histograms_FM2.csv")

# Drop non-feature columns
features = [f for f in df.columns if f not in ["label", "weight", "process"]]
X = df[features]

# Compute correlation matrix
corr = X.corr()

# Rename columns and index with LaTeX labels
corr_latex = corr.rename(index=feature_labels, columns=feature_labels)

# Plot heatmap
plt.figure(figsize=(14, 12))
sns.heatmap(corr_latex, cmap="coolwarm", square=True,
            cbar_kws={"label": "Correlation"}, annot=False, fmt=".2f",
            xticklabels=True, yticklabels=True)

plt.xticks(rotation=45, ha="right")
plt.yticks(rotation=0)
plt.title("Feature Correlation Matrix (LaTeX labels)", fontsize=16)
plt.tight_layout()
plt.savefig("feature_correlation_matrix_latex.pdf")
plt.show()
