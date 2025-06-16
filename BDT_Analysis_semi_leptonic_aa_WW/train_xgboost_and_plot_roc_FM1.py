

## ================================================================================
##        Hamzeh Khanpour  --- June 2025
## ================================================================================



import pandas as pd
import xgboost as xgb
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_curve, auc
from sklearn.metrics import roc_curve, roc_auc_score

import mplhep as hep
#hep.style.use("CMS")


# Load dataset
df = pd.read_csv("ml_input_from_histograms_FM1.csv")

# Extract label, weight, process before splitting
y = df["label"]
weights = df["weight"]
process = df["process"]

# Remove non-feature columns
X = df.drop(columns=["label", "weight", "process"])

# Split including process
X_train, X_test, y_train, y_test, w_train, w_test, proc_train, proc_test = train_test_split(
    X, y, weights, process, test_size=0.25, stratify=y, random_state=42
)

# Train XGBoost with weights
model = xgb.XGBClassifier(
    n_estimators=200,
    max_depth=4,
    learning_rate=0.05,
    use_label_encoder=False,
    eval_metric="logloss"
)
model.fit(X_train, y_train, sample_weight=w_train)

# Predict probabilities
y_scores = model.predict_proba(X_test)[:, 1]

# Save for limit setting
df_out = pd.DataFrame({
    "bdt_score": y_scores,
    "label": y_test.values,
    "weight": w_test.values,
    "process": proc_test.values
})
df_out.to_csv("ml_with_bdt_scores_FM1.csv", index=False)
print("✅ Saved: ml_with_bdt_scores_FM1.csv with BDT scores, labels, weights, and processes.")


fpr, tpr, _ = roc_curve(y_test, y_scores)
roc_auc = roc_auc_score(y_test, y_scores)


plt.figure(figsize=(8, 6))


plt.plot(fpr, tpr, label=f"AUC = {roc_auc:.4f}", linewidth=2)
plt.plot([0, 1], [0, 1], 'k--', alpha=0.5)
plt.xlabel("False Positive Rate")
plt.ylabel("True Positive Rate")
plt.title("ROC Curve (Signal vs Background)")
plt.legend(loc="lower right")
plt.grid(True, linestyle="--", alpha=0.5)
plt.tight_layout()
plt.savefig("roc_curve_xgboost_FM1.pdf")
plt.show()



# Feature importances
plt.figure(figsize=(8, 6))
xgb.plot_importance(model, importance_type="gain", show_values=False)
plt.title("Feature Importances (Gain)")
plt.tight_layout()
plt.savefig("feature_importances_xgboost_FM1.pdf")
plt.show()



