
import numpy as np
import pandas as pd
import xgboost as xgb
from sklearn.model_selection import train_test_split

# Load prepared CSV data
df = pd.read_csv("ml_input_from_histograms.csv")

# Features, labels, weights, process tags
X = df.drop(columns=["label", "weight", "process"])
y = df["label"]
w = df["weight"]
process = df["process"]

# Train/test split with stratification
X_train, X_test, y_train, y_test, w_train, w_test, proc_train, proc_test = train_test_split(
    X, y, w, process, test_size=0.3, stratify=y, random_state=42
)

# Load or train model
model = xgb.XGBClassifier(n_estimators=200, max_depth=4, learning_rate=0.05, use_label_encoder=False, eval_metric='logloss')
model.fit(X_train, y_train, sample_weight=w_train)

# Apply BDT model
y_scores = model.predict_proba(X_test)[:, 1]

# Define BDT score cut
cut = 0.232
pass_cut = y_scores > cut

# Signal and background masks
signal_mask = y_test == 1
background_mask = y_test == 0

# Signal and background efficiencies
signal_eff = np.sum(w_test[signal_mask & pass_cut]) / np.sum(w_test[signal_mask])
background_eff = np.sum(w_test[background_mask & pass_cut]) / np.sum(w_test[background_mask])

# Integrated luminosity
luminosity_pb = 100.0 * 1000.0  # [pb^-1]

# Expected yields at 100 fb^-1 from BDT test sample
expected_signal_yield = np.sum(w_test[signal_mask & pass_cut])
expected_background_yield = np.sum(w_test[background_mask & pass_cut])

# Breakdown by background process (from test set)
background_df = pd.DataFrame({
    "process": proc_test[background_mask & pass_cut],
    "weight": w_test[background_mask & pass_cut]
})
bkg_yields_by_process = background_df.groupby("process")["weight"].sum().sort_values(ascending=False)

# Cross-section × efficiency method
cross_section_eff = {
    "aa_ww": [0.0099465, 0.0315],
    "aa_ttbar": [4.824851e-05, 0.0074],
    "aa_tautau": [2.51510000, 0.0000],
    "inclusive_ttbar": [0.0065764, 0.0019],
    "w_production": [1.910288, 0.0001],
    "z_production": [0.24064758729900002, 0.0001],
    "wwj": [0.016080595320336195, 0.0020],
    "zzj": [6.694889944457796e-05, 0.0029],
    "wzj": [0.0023785292894910495, 0.0018],
}

expected_bkg_physics = {k: sigma * eff * luminosity_pb for k, (sigma, eff) in cross_section_eff.items()}
expected_bkg_total = sum(expected_bkg_physics.values())

# Output results
print(f"\n📊 Post-BDT Cut Summary (Score > {cut})")
print(f"Signal Efficiency: {signal_eff:.4f}")
print(f"Background Efficiency (test set weights): {background_eff:.4f}")
print(f"Expected Signal Yield (100 fb⁻¹, weighted test set): {expected_signal_yield:.2f}")
print(f"Expected Background Yield (100 fb⁻¹, weighted test set): {expected_background_yield:.2f}")

print(f"\n✅ Physics-Based Expected Background (σ × ε × L): {expected_bkg_total:.2f} events")
print("\n🧾 Background Composition (test set weighted):")
print(bkg_yields_by_process.to_string())
