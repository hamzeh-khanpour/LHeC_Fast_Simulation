import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_curve, auc, accuracy_score
from sklearn.preprocessing import StandardScaler

# Classifiers
from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import RandomForestClassifier, AdaBoostClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.svm import SVC
from sklearn.neural_network import MLPClassifier
import xgboost as xgb

# Load dataset
df = pd.read_csv("ml_input_from_histograms_FM2.csv")

drop_cols = ["label", "weight"]
if "process" in df.columns:
    drop_cols.append("process")

X = df.drop(columns=drop_cols)
y = df["label"]
weights = df["weight"]


# Train/test split
X_train, X_test, y_train, y_test, w_train, w_test = train_test_split(
    X, y, weights, stratify=y, test_size=0.3, random_state=42
)

# Normalize
scaler = StandardScaler()
X_train_scaled = scaler.fit_transform(X_train)
X_test_scaled = scaler.transform(X_test)

# Models
models = {
    "DecisionTree": DecisionTreeClassifier(),
    "RandomForest": RandomForestClassifier(n_estimators=100),
    "AdaBoost": AdaBoostClassifier(n_estimators=200),
    "KNN": KNeighborsClassifier(n_neighbors=5),
    "SVM": SVC(probability=True, kernel="rbf"),
    "MLP": MLPClassifier(hidden_layer_sizes=(50,), max_iter=1000),
    "TMVABDT": xgb.XGBClassifier(
        n_estimators=850, learning_rate=0.5, max_depth=3, subsample=0.5,
        use_label_encoder=False, eval_metric="logloss"
    )
}

# Classifiers that do NOT support sample_weight
no_weight_support = ["KNN", "SVM", "MLP"]

# ROC curve plot
plt.figure(figsize=(10, 8))

for name, model in models.items():
    print(f"🔧 Training {name}...")

    if name in no_weight_support:
        model.fit(X_train_scaled, y_train)
    else:
        model.fit(X_train_scaled, y_train, sample_weight=w_train)

    if hasattr(model, "predict_proba"):
        y_scores = model.predict_proba(X_test_scaled)[:, 1]
    else:
        y_scores = model.decision_function(X_test_scaled)

    # ROC and AUC with test weights
    fpr, tpr, _ = roc_curve(y_test, y_scores, sample_weight=w_test)
    auc_val = auc(fpr, tpr)

    # Accuracy (unweighted)
    y_pred = model.predict(X_test_scaled)
    acc = accuracy_score(y_test, y_pred)

    plt.plot(fpr, tpr, label=f"{name} (AUC = {auc_val:.3f})")
    print(f"✅ {name} Accuracy: {acc:.4f}")

# Finalize ROC plot
plt.plot([0, 1], [0, 1], 'k--', alpha=0.5)
plt.xlabel("False Positive Rate")
plt.ylabel("True Positive Rate")
plt.title("ROC Curve Comparison: Signal vs Background")
plt.legend()
plt.grid(True, linestyle="--", alpha=0.5)
plt.tight_layout()
plt.savefig("ml_comparison_roc_curves_FM2.pdf")
plt.show()
