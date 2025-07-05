## ================================================================================
##     Unnormalized Distributions with S/B Ratio Panel After BDT Cut
##     Hamzeh Khanpour — Updated by ChatGPT (July 2025)
## ================================================================================

import pandas as pd
import matplotlib.pyplot as plt
import xgboost as xgb
import numpy as np
import matplotlib.gridspec as gridspec

# Parameters
BDT_THRESHOLD = 0.253  # Best BDT cut from optimization

# Load dataset
df = pd.read_csv("ml_input_from_histograms_FM2.csv")

# Drop non-feature columns
drop_cols = ["label", "weight"]
if "process" in df.columns:
    drop_cols.append("process")

X = df.drop(columns=drop_cols)
y = df["label"]
weights = df["weight"]

# Train BDT
model = xgb.XGBClassifier(
    n_estimators=200,
    max_depth=4,
    learning_rate=0.05,
    use_label_encoder=False,
    eval_metric="logloss"
)
model.fit(X, y, sample_weight=weights)

# Apply BDT cut
df["bdt_score"] = model.predict_proba(X)[:, 1]
df_selected = df[df["bdt_score"] > BDT_THRESHOLD]
print(f"💘 Events surviving BDT > {BDT_THRESHOLD}: {len(df_selected)}")

# Features to analyze
features_to_plot = [
    "lepton_pt", "missing_et", "leading_jet_pt",
    "m_w_leptonic", "m_w_hadronic", "m_lvjj", "ht_total", "mt_w_leptonic"
]

# Plot unnormalized distributions and S/B
for feature in features_to_plot:
    fig = plt.figure(figsize=(8, 7))
    gs = gridspec.GridSpec(2, 1, height_ratios=[3, 1], hspace=0.05)

    # Main Plot
    ax_main = plt.subplot(gs[0])
    bins = np.linspace(
        df_selected[feature].min(),
        df_selected[feature].max(),
        50
    )

    s_vals, _ = np.histogram(
        df_selected[df_selected["label"] == 1][feature],
        bins=bins,
        weights=df_selected[df_selected["label"] == 1]["weight"]
    )
    b_vals, _ = np.histogram(
        df_selected[df_selected["label"] == 0][feature],
        bins=bins,
        weights=df_selected[df_selected["label"] == 0]["weight"]
    )

    bin_centers = 0.5 * (bins[1:] + bins[:-1])

    ax_main.hist(df[df["label"] == 1][feature], bins=bins,
                 weights=df[df["label"] == 1]["weight"], alpha=0.3,
                 label="Signal", color="red")

    ax_main.hist(df[df["label"] == 0][feature], bins=bins,
                 weights=df[df["label"] == 0]["weight"], alpha=0.3,
                 label="Background", color="blue")

    ax_main.step(bin_centers, s_vals, where='mid', color='darkred',
                 label="Signal (BDT cut)", linewidth=2)
    ax_main.step(bin_centers, b_vals, where='mid', color='navy',
                 label="Background (BDT cut)", linewidth=2)

    ax_main.set_ylabel("Events")
    ax_main.set_title(f"{feature} Distribution")
    ax_main.legend()
    ax_main.set_yscale("log")
    ax_main.grid(True, linestyle="--", alpha=0.5)

    # Ratio Plot
    ax_ratio = plt.subplot(gs[1], sharex=ax_main)
    with np.errstate(divide='ignore', invalid='ignore'):
        ratio = np.divide(s_vals, b_vals, out=np.zeros_like(s_vals), where=b_vals > 0)
    ax_ratio.step(bin_centers, ratio, where='mid', color='black', linewidth=2)
    ax_ratio.set_ylabel("S/B", fontsize=11)
    ax_ratio.set_xlabel(feature)
    ax_ratio.grid(True, linestyle="--", alpha=0.5)
    ax_ratio.set_ylim(0, np.max(ratio) * 1.2 if np.max(ratio) > 0 else 1)

    plt.tight_layout()
    plt.savefig(f"unnormalized_{feature}_with_SB_ratio_bdt_cut_FM2.pdf")
    plt.show()

print("✅ All unnormalized plots with S/B ratio panels saved.")


