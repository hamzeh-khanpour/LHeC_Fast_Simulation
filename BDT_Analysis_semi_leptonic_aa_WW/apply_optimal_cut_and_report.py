import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import xgboost as xgb

# Load dataset
df = pd.read_csv("ml_input_from_histograms.csv")
X = df.drop(columns=["label"])
y = df["label"]

# Train model
model = xgb.XGBClassifier(
    n_estimators=200,
    max_depth=4,
    learning_rate=0.05,
    use_label_encoder=False,
    eval_metric="logloss"
)
model.fit(X, y)

# Predict BDT scores
df["bdt_score"] = model.predict_proba(X)[:, 1]

# Find optimal BDT threshold
thresholds = np.linspace(0.0, 1.0, 100)
best_threshold = 0.0
best_significance = 0.0

for t in thresholds:
    selected = df[df["bdt_score"] > t]
    s = len(selected[selected["label"] == 1])
    b = len(selected[selected["label"] == 0])
    if s + b > 0:
        Z = s / np.sqrt(s + b)
        if Z > best_significance:
            best_significance = Z
            best_threshold = t

# Apply best cut
df_cut = df[df["bdt_score"] > best_threshold]

# Report
s_total = len(df[df["label"] == 1])
s_after = len(df_cut[df_cut["label"] == 1])
b_after = len(df_cut[df_cut["label"] == 0])

efficiency = s_after / s_total if s_total > 0 else 0.0

print(f"💎 Optimal BDT cut: {best_threshold:.3f}")
print(f"🧬 Signal efficiency: {efficiency:.4f} ({s_after}/{s_total})")
print(f"🎭 Background after cut: {b_after} events")

# Save to CSV
df_cut[df_cut["label"] == 1].to_csv("signal_after_cut.csv", index=False)
df_cut[df_cut["label"] == 0].to_csv("background_after_cut.csv", index=False)

# Optional: plot score distribution with cut
plt.figure(figsize=(8, 6))
plt.hist(df[df["label"] == 1]["bdt_score"], bins=50, alpha=0.5, label="Signal", density=True)
plt.hist(df[df["label"] == 0]["bdt_score"], bins=50, alpha=0.5, label="Background", density=True)
plt.axvline(best_threshold, color='red', linestyle='--', label=f"Cut = {best_threshold:.3f}")
plt.xlabel("BDT Score")
plt.ylabel("Normalized Entries")
plt.title("BDT Score with Optimal Cut")
plt.legend()
plt.grid(True, linestyle="--", alpha=0.5)
plt.tight_layout()
plt.savefig("bdt_cut_applied.pdf")
plt.show()
