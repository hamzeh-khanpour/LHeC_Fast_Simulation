import pandas as pd
import xgboost as xgb
import shap
import matplotlib.pyplot as plt

# Load data
df = pd.read_csv("ml_input_from_histograms_FM2.csv")

drop_cols = ["label"]
if "process" in df.columns:
    drop_cols.append("process")
X = df.drop(columns=drop_cols)
y = df["label"]

# Train model
#model = xgb.XGBClassifier(
    #n_estimators=200,
    #max_depth=4,
    #learning_rate=0.05,
    #use_label_encoder=False,
    #eval_metric="logloss"
#)

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


model.fit(X, y)

# Explain with SHAP
explainer = shap.Explainer(model)
shap_values = explainer(X)

# 1️⃣ Summary plot: Feature importance
plt.figure()
shap.summary_plot(shap_values, X, show=False)
plt.tight_layout()
plt.savefig("shap_summary_plot_FM2.pdf")
plt.close()

# 2️⃣ Dependence plot
shap.dependence_plot(
    "jet_centrality", shap_values.values, X,
    interaction_index="leading_jet_eta", show=False
)
plt.title("SHAP: jet_centrality vs leading_jet_eta")
plt.tight_layout()
plt.savefig("shap_jet_centrality_vs_leading_eta_FM2.pdf")
plt.close()

print("🌿 SHAP analysis complete. Files saved:")
print("  ➤ shap_summary_plot.pdf")
print("  ➤ shap_jet_centrality_vs_leading_eta.pdf")
