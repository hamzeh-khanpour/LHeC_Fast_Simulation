## ================================================================================
##     Unnormalized Distributions Before/After BDT Cut
##     Hamzeh Khanpour — Updated by ChatGPT
## ================================================================================

import pandas as pd
import matplotlib.pyplot as plt
import xgboost as xgb

# Parameters
BDT_THRESHOLD = 0.253  # This should ideally come from optimization

# Load dataset
df = pd.read_csv("ml_input_from_histograms_FM2.csv")

# Drop non-feature columns
drop_cols = ["label", "weight"]
if "process" in df.columns:
    drop_cols.append("process")

X = df.drop(columns=drop_cols)
y = df["label"]
weights = df["weight"]

# Retrain XGBoost with weights (for consistency)
model = xgb.XGBClassifier(
    n_estimators=200,
    max_depth=4,
    learning_rate=0.05,
    use_label_encoder=False,
    eval_metric="logloss"
)
model.fit(X, y, sample_weight=weights)

# Predict and apply BDT cut
df["bdt_score"] = model.predict_proba(X)[:, 1]
df_selected = df[df["bdt_score"] > BDT_THRESHOLD]
print(f"💘 Events surviving BDT > {BDT_THRESHOLD}: {len(df_selected)}")

# Features to plot
features_to_plot = [
    "lepton_pt", "missing_et", "leading_jet_pt", "m_w_leptonic",
    "m_w_hadronic", "m_lvjj", "ht_total", "mt_w_leptonic"
]

# Plot raw (unnormalized) yields before/after cut
for feature in features_to_plot:
    plt.figure(figsize=(8, 6))

    # Raw yields before BDT
    plt.hist(df[df["label"] == 1][feature], bins=50,
             weights=df[df["label"] == 1]["weight"],
             alpha=0.4, label="Signal", color="red")

    plt.hist(df[df["label"] == 0][feature], bins=50,
             weights=df[df["label"] == 0]["weight"],
             alpha=0.4, label="Background", color="blue")

    # Raw yields after BDT
    plt.hist(df_selected[df_selected["label"] == 1][feature], bins=50,
             weights=df_selected[df_selected["label"] == 1]["weight"],
             alpha=0.6, label="Signal (BDT cut)", histtype='stepfilled', edgecolor='darkred', color="red")

    plt.hist(df_selected[df_selected["label"] == 0][feature], bins=50,
             weights=df_selected[df_selected["label"] == 0]["weight"],
             alpha=0.6, label="Background (BDT cut)", histtype='stepfilled', edgecolor='navy', color="blue")

    plt.title(f"{feature} Distribution (Unnormalized Yields)")
    plt.xlabel(feature)
    plt.ylabel("Events")
    plt.legend()
    plt.grid(True, linestyle="--", alpha=0.5)
    plt.tight_layout()
    plt.savefig(f"unnormalized_{feature}_bdt_cut_FM2.pdf")
    plt.yscale("log")
    plt.show()

print("✅ All unnormalized plots saved with BDT selection overlay.")
