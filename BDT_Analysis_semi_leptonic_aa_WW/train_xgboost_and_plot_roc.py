import pandas as pd
import xgboost as xgb
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_curve, auc

# Load dataset
df = pd.read_csv("ml_input_from_histograms.csv")
X = df.drop(columns=["label", "weight", "process"])

y = df["label"]
weights = df["weight"]

# Split
X_train, X_test, y_train, y_test, w_train, w_test = train_test_split(
    X, y, weights, test_size=0.25, stratify=y, random_state=42
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

# ROC curve
fpr, tpr, _ = roc_curve(y_test, y_scores)
roc_auc = auc(fpr, tpr)

plt.figure(figsize=(8, 6))
plt.plot(fpr, tpr, label=f"AUC = {roc_auc:.4f}", linewidth=2)
plt.plot([0, 1], [0, 1], 'k--', alpha=0.5)
plt.xlabel("False Positive Rate")
plt.ylabel("True Positive Rate")
plt.title("ROC Curve (Signal vs Background)")
plt.legend(loc="lower right")
plt.grid(True, linestyle="--", alpha=0.5)
plt.tight_layout()
plt.savefig("roc_curve_xgboost.pdf")
plt.show()

# Feature importances
plt.figure(figsize=(10, 6))
xgb.plot_importance(model, importance_type="gain", show_values=False)
plt.title("Feature Importances (Gain)")
plt.tight_layout()
plt.savefig("feature_importances_xgboost.pdf")
plt.show()

# BDT score distribution
plt.figure(figsize=(8, 6))
plt.hist(y_scores[y_test == 1], bins=50, alpha=0.6, label="Signal", density=True)
plt.hist(y_scores[y_test == 0], bins=50, alpha=0.6, label="Background", density=True)
plt.xlabel("XGBoost BDT Score")
plt.ylabel("Normalized Events")
plt.title("BDT Score Distribution")
plt.legend()
plt.grid(True, linestyle="--", alpha=0.5)
plt.tight_layout()
plt.savefig("bdt_score_distribution.pdf")
plt.show()
