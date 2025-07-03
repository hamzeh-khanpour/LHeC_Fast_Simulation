import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

# Load the pre-cut ML input
df = pd.read_csv("ml_input_from_histograms_FM2.csv")

# Drop any non-numerical metadata columns
drop_cols = ["label", "weight"]
if "process" in df.columns:
    drop_cols.append("process")

X = df.drop(columns=drop_cols)
weights = df["weight"]
labels = df["label"]


# 1️⃣ Full Correlation Heatmap
plt.figure(figsize=(12, 10))
sns.heatmap(X.corr(), cmap="coolwarm", annot=True, fmt=".2f", square=True)
plt.title("🔍 Full Correlation Heatmap")
plt.tight_layout()
plt.savefig("diagnostic_full_corr_FM2.pdf")
plt.close()

# 2️⃣ Jet-Centrality Focused Scatter Plots
plt.figure(figsize=(6, 5))
sns.scatterplot(data=df, x="leading_jet_eta", y="jet_centrality", hue="label", alpha=0.5, palette=["blue", "red"])
plt.title("Jet Centrality vs Leading Jet Eta")
plt.tight_layout()
plt.savefig("scatter_jet_centrality_vs_leading_eta_FM2.pdf")
plt.close()

plt.figure(figsize=(6, 5))
sns.scatterplot(data=df, x="subleading_jet_eta", y="jet_centrality", hue="label", alpha=0.5, palette=["blue", "red"])
plt.title("Jet Centrality vs Subleading Jet Eta")
plt.tight_layout()
plt.savefig("scatter_jet_centrality_vs_subleading_eta_FM2.pdf")
plt.close()

# 3️⃣ Principal Component Analysis (PCA)
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)

pca = PCA()
X_pca = pca.fit_transform(X_scaled)
explained_var = pca.explained_variance_ratio_

plt.figure(figsize=(8, 6))
plt.plot(range(1, len(explained_var) + 1), explained_var, marker='o')
plt.xlabel("Principal Component")
plt.ylabel("Variance Explained")
plt.title("PCA: Variance Explained by Components")
plt.grid(True, linestyle='--', alpha=0.5)
plt.tight_layout()
plt.savefig("pca_variance_explained_FM2.pdf")
plt.close()

# 4️⃣ Print strongest correlated pairs (above 0.8)
print("\n🔍 Strongly Correlated Feature Pairs (|corr| > 0.8):")
corr_matrix = X.corr().abs()
for i in range(len(corr_matrix.columns)):
    for j in range(i + 1, len(corr_matrix.columns)):
        corr_val = corr_matrix.iloc[i, j]
        if corr_val > 0.8:
            f1 = corr_matrix.columns[i]
            f2 = corr_matrix.columns[j]
            print(f"{f1} ↔ {f2} = {corr_val:.3f}")
