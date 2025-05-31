import pandas as pd
import matplotlib.pyplot as plt
import xgboost as xgb

# Parameters
BDT_THRESHOLD = 0.374  # ← Ideally retrieved from optimize_bdt_cut.py

# Load dataset
df = pd.read_csv("ml_input_from_histograms.csv")

# Drop any categorical columns if they exist
drop_cols = ["label", "weight"]
if "process" in df.columns:
    drop_cols.append("process")

X = df.drop(columns=drop_cols)
y = df["label"]
weights = df["weight"]

# Retrain the model with weights (to ensure consistency)
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
    "lepton_pt", "missing_et", "jet_centrality", "m_w_leptonic"
]

# Plot before/after cut using event weights
for feature in features_to_plot:
    plt.figure(figsize=(8, 6))

    # Raw distributions
    plt.hist(df[df["label"] == 1][feature], bins=50, weights=df[df["label"] == 1]["weight"],
             density=True, alpha=0.4, label="Signal (raw)")
    plt.hist(df[df["label"] == 0][feature], bins=50, weights=df[df["label"] == 0]["weight"],
             density=True, alpha=0.4, label="Background (raw)")

    # After BDT cut
    plt.hist(df_selected[df_selected["label"] == 1][feature], bins=50,
             weights=df_selected[df_selected["label"] == 1]["weight"],
             density=True, alpha=0.6, label="Signal (BDT cut)", histtype='stepfilled', edgecolor='red')
    plt.hist(df_selected[df_selected["label"] == 0][feature], bins=50,
             weights=df_selected[df_selected["label"] == 0]["weight"],
             density=True, alpha=0.6, label="Background (BDT cut)", histtype='stepfilled', edgecolor='blue')

    plt.title(f"{feature} Distribution Before/After BDT Cut")
    plt.xlabel(feature)
    plt.ylabel("Normalized Entries")
    plt.legend()
    plt.grid(True, linestyle="--", alpha=0.5)
    plt.tight_layout()
    plt.savefig(f"{feature}_bdt_cut_comparison.pdf")

print("✅ All plots saved with BDT selection overlay.")
