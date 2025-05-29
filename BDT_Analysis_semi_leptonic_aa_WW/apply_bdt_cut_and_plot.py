import pandas as pd
import matplotlib.pyplot as plt
import xgboost as xgb

# Parameters
BDT_THRESHOLD = 0.384  # Feel free to optimize later

# Load dataset
df = pd.read_csv("ml_input_from_histograms.csv")
X = df.drop(columns=["label"])
y = df["label"]

# Load or re-train the model
model = xgb.XGBClassifier(
    n_estimators=200,
    max_depth=4,
    learning_rate=0.05,
    use_label_encoder=False,
    eval_metric="logloss"
)
model.fit(X, y)

# Compute BDT scores
df["bdt_score"] = model.predict_proba(X)[:, 1]

# Apply BDT score cut
df_selected = df[df["bdt_score"] > BDT_THRESHOLD]
print(f"💘 Events surviving BDT > {BDT_THRESHOLD}: {len(df_selected)}")

# Define features to replot
features_to_plot = [
    "lepton_pt", "missing_et", "jet_centrality", "m_w_leptonic"
]

# Plot comparison: original vs selected
for feature in features_to_plot:
    plt.figure(figsize=(8, 6))

    # Original distributions
    plt.hist(df[df["label"] == 1][feature], bins=50, density=True, alpha=0.4, label="Signal (raw)")
    plt.hist(df[df["label"] == 0][feature], bins=50, density=True, alpha=0.4, label="Background (raw)")

    # After BDT cut
    plt.hist(df_selected[df_selected["label"] == 1][feature], bins=50, density=True,
             alpha=0.6, label="Signal (BDT cut)", histtype='stepfilled', edgecolor='red')
    plt.hist(df_selected[df_selected["label"] == 0][feature], bins=50, density=True,
             alpha=0.6, label="Background (BDT cut)", histtype='stepfilled', edgecolor='blue')

    plt.title(f"{feature} Distribution Before/After BDT Cut")
    plt.xlabel(feature)
    plt.ylabel("Normalized Entries")
    plt.legend()
    plt.grid(True, linestyle="--", alpha=0.5)
    plt.tight_layout()
    plt.savefig(f"{feature}_bdt_cut_comparison.pdf")

print("✅ All plots saved with BDT selection overlay.")
