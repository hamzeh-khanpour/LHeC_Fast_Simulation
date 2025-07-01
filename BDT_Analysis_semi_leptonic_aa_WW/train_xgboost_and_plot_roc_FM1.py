

## ================================================================================
##        Hamzeh Khanpour  --- June 2025
## ================================================================================


import pandas as pd
import xgboost as xgb
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_curve, auc
from sklearn.metrics import roc_curve, roc_auc_score
import seaborn as sns
import numpy as np


from sklearn.metrics import (
    roc_curve, roc_auc_score,
    precision_recall_curve, average_precision_score
)



import mplhep as hep

# Use CMS plot style
hep.style.use("CMS")


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
    n_estimators=500,
    max_depth=5,
    learning_rate=0.02,
    subsample=0.8,
    colsample_bytree=0.8,
    min_child_weight=2,
    use_label_encoder=False,
    eval_metric="auc"
)


## Train XGBoost with weights
#model = xgb.XGBClassifier(
    #n_estimators=200,
    #max_depth=4,
    #learning_rate=0.05,
    #use_label_encoder=False,
    #eval_metric="logloss"
#)



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



# ------------------------------
# Predict BDT scores
# ------------------------------
y_scores_test = model.predict_proba(X_test)[:, 1]
y_scores_train = model.predict_proba(X_train)[:, 1]



plt.figure(figsize=(11, 12))
plt.subplots_adjust(left=0.15, right=0.95, bottom=0.12, top=0.95)



# ------------------------------
# 1. BDT Score Distribution (Signal vs Background)
# ------------------------------

plt.hist(y_scores_test[y_test == 1], bins=50, weights=w_test[y_test == 1],
         alpha=0.6, label="Signal", density=True)
plt.hist(y_scores_test[y_test == 0], bins=50, weights=w_test[y_test == 0],
         alpha=0.6, label="Background", density=True)
plt.xlabel("BDT Score")
plt.ylabel("Normalized Events")
plt.title("BDT Score Distribution (Test Set)")
plt.legend()
plt.grid(True, linestyle="--", alpha=0.5)
plt.tight_layout()
plt.savefig("bdt_score_distribution_test.pdf")
plt.show()



# ------------------------------
# 2. ROC Curve
# ------------------------------
fpr, tpr, _ = roc_curve(y_test, y_scores_test)
roc_auc = roc_auc_score(y_test, y_scores_test)

plt.plot(fpr, tpr, label=f"AUC = {roc_auc:.4f}", linewidth=2)
plt.plot([0, 1], [0, 1], "k--", alpha=0.5)
plt.xlabel("False Positive Rate")
plt.ylabel("True Positive Rate")
plt.title("ROC Curve")
plt.legend(loc="lower right")
plt.grid(True, linestyle="--", alpha=0.5)
plt.tight_layout()
plt.savefig("roc_curve.pdf")
plt.show()



# ------------------------------
# 3. Precision-Recall Curve
# ------------------------------
precision, recall, _ = precision_recall_curve(y_test, y_scores_test)
ap = average_precision_score(y_test, y_scores_test)

plt.plot(recall, precision, label=f"AP = {ap:.4f}", linewidth=2)
plt.xlabel("Recall")
plt.ylabel("Precision")
plt.title("Precision-Recall Curve")
plt.legend()
plt.grid(True, linestyle="--", alpha=0.5)
plt.tight_layout()
plt.savefig("precision_recall_curve.pdf")
plt.show()



# ------------------------------
# 4. Feature Importance (Gain)
# ------------------------------
xgb.plot_importance(model, importance_type="gain", show_values=False)
plt.title("Feature Importances (Gain)")
plt.tight_layout()
plt.savefig("feature_importances_gain.pdf")
plt.show()




# ------------------------------
# 5. Feature Correlation Matrix
# ------------------------------
corr = X.corr()
sns.heatmap(corr, cmap="coolwarm", square=True, cbar_kws={'label': 'Correlation'})
plt.title("Correlation Matrix of Input Features")
plt.tight_layout()
plt.savefig("feature_correlation_matrix.pdf")
plt.show()



# ------------------------------
# 6. BDT Score by Background Process
# ------------------------------
for proc_name in np.unique(proc_test):
    mask = (proc_test == proc_name) & (y_test == 0)
    plt.hist(y_scores_test[mask], bins=50, weights=w_test[mask],
             label=proc_name, histtype="step", density=True)
plt.hist(y_scores_test[y_test == 1], bins=50, weights=w_test[y_test == 1],
         label="Signal", histtype="step", linewidth=2, color="black")
plt.xlabel("BDT Score")
plt.ylabel("Normalized Events")
plt.title("BDT Score by Background Process")
plt.legend()
plt.grid(True, linestyle="--", alpha=0.5)
plt.tight_layout()
plt.savefig("bdt_score_by_process.pdf")
plt.show()

print("✅ All BDT diagnostic plots saved successfully.")





# Define range of BDT score thresholds
bdt_cut_bins = np.linspace(0, 1, 200)
significance_values = []
bdt_thresholds = []

# Loop over BDT cut values
for cut in bdt_cut_bins:
    selected = df_out[df_out["bdt_score"] >= cut]

    S = selected[selected["label"] == 1]["weight"].sum()
    B = selected[selected["label"] == 0]["weight"].sum()

    if (S + B) > 0:
        significance = S / np.sqrt(S + B)
    else:
        significance = 0

    bdt_thresholds.append(cut)
    significance_values.append(significance)

# Get max significance point
max_significance = max(significance_values)
optimal_cut = bdt_thresholds[np.argmax(significance_values)]

# Plot
plt.figure(figsize=(10, 6))
plt.plot(bdt_thresholds, significance_values, label="S / √(S + B)", linewidth=2)
plt.axvline(optimal_cut, color='r', linestyle='--', label=f"Optimal Cut = {optimal_cut:.3f}")
plt.xlabel("BDT Score Threshold")
plt.ylabel("Significance (S / √(S + B))")
plt.title("Signal Significance vs BDT Cut")
plt.grid(True, linestyle="--", alpha=0.5)
plt.legend()
plt.tight_layout()
plt.savefig("significance_vs_bdt_cut.pdf")
plt.show()





