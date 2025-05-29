import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import xgboost as xgb

# Load dataset
df = pd.read_csv("ml_input_from_histograms.csv")
X = df.drop(columns=["label"])
y = df["label"]

# Train XGBoost model
model = xgb.XGBClassifier(
    n_estimators=200,
    max_depth=4,
    learning_rate=0.05,
    use_label_encoder=False,
    eval_metric="logloss"
)
model.fit(X, y)

# Get BDT scores
df["bdt_score"] = model.predict_proba(X)[:, 1]

# Scan thresholds
thresholds = np.linspace(0.0, 1.0, 100)
significances = []

for t in thresholds:
    selected = df[df["bdt_score"] > t]
    s = len(selected[selected["label"] == 1])
    b = len(selected[selected["label"] == 0])
    if s + b > 0:
        Z = s / np.sqrt(s + b)
    else:
        Z = 0
    significances.append(Z)

# Find optimal threshold
best_index = np.argmax(significances)
best_threshold = thresholds[best_index]
best_significance = significances[best_index]

# Plot
plt.figure(figsize=(8, 6))
plt.plot(thresholds, significances, label="Z = S / sqrt(S + B)")
plt.axvline(best_threshold, color='red', linestyle='--',
            label=f"Optimal cut = {best_threshold:.3f}")
plt.title("BDT Cut Optimization")
plt.xlabel("BDT Score Threshold")
plt.ylabel("Significance (Z)")
plt.grid(True, linestyle="--", alpha=0.5)
plt.legend()
plt.tight_layout()
plt.savefig("optimal_bdt_cut.pdf")
plt.show()

print(f"💘 Best threshold: {best_threshold:.3f} with significance Z = {best_significance:.3f}")
