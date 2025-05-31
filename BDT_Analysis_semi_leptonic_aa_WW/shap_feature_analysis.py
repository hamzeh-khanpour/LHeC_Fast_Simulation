import pandas as pd
import xgboost as xgb
import shap
import matplotlib.pyplot as plt

# Load data
df = pd.read_csv("ml_input_from_histograms.csv")

drop_cols = ["label"]
if "process" in df.columns:
    drop_cols.append("process")
X = df.drop(columns=drop_cols)
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

# Explain with SHAP
explainer = shap.Explainer(model)
shap_values = explainer(X)

# 1️⃣ Summary plot: Feature importance
plt.figure()
shap.summary_plot(shap_values, X, show=False)
plt.tight_layout()
plt.savefig("shap_summary_plot.pdf")
plt.close()

# 2️⃣ Dependence plot
shap.dependence_plot(
    "jet_centrality", shap_values.values, X,
    interaction_index="leading_jet_eta", show=False
)
plt.title("SHAP: jet_centrality vs leading_jet_eta")
plt.tight_layout()
plt.savefig("shap_jet_centrality_vs_leading_eta.pdf")
plt.close()

print("🌿 SHAP analysis complete. Files saved:")
print("  ➤ shap_summary_plot.pdf")
print("  ➤ shap_jet_centrality_vs_leading_eta.pdf")
