
# 📊 Anomalous Quartic Gauge Couplings (aQGC) Sensitivity at the LHeC

This repository contains a full pipeline for estimating the sensitivity of the LHeC to the anomalous quartic gauge coupling parameter $f_{M2}/\Lambda^4$ via shape-based profile likelihood scans. The analysis is based on semi-leptonic $W^+W^-$ production from $\gamma\gamma$ fusion, including realistic event simulation, machine learning-based signal-background discrimination, and statistical limit setting.

---

## 📁 Project Structure

```
├── ExtractData_Train_PlotROC.py           # Extracts and saves histogram-based ML features
├── train_xgboost_and_plot_roc.py          # Trains a BDT and computes ROC + saves scores
├── limit_scan_fm2_profile_likelihood_v5.py  # Full likelihood-based limit calculation
├── output_histograms.root                 # ROOT file with preselection histograms
├── ml_input_from_histograms.csv           # Raw ML inputs from histograms
├── ml_with_bdt_scores.csv                 # ML inputs with predicted BDT scores
```

---

## 🧠 ML Pipeline

### 1. Data Preparation

Run:

```bash
python3 ExtractData_Train_PlotROC.py
```

- Extracts samples from histogram templates in `output_histograms.root`
- Samples kinematic variables from:
  - Leptons
  - Jets
  - MET
- Saves as `ml_input_from_histograms.csv` with `label`, `weight`, and `process`.

### 2. Training the BDT Classifier

Run:

```bash
python3 train_xgboost_and_plot_roc.py
```

- Trains XGBoost BDT
- Computes AUC & ROC curve
- Appends `bdt_score` and `process` columns
- Saves as `ml_with_bdt_scores.csv`

---

## 📐 Profile Likelihood Limit Setting

Run:

```bash
python3 limit_scan_fm2_profile_likelihood_v5.py
```

### Features:

- Calculates optimal BDT cut based on AMS.
- Computes:
  - ML efficiency
  - Preselection efficiency
  - Total efficiency = ML × Preselection
- Constructs Asimov dataset (background only).
- Defines likelihood:

  $$
  \mathcal{L}(f_{M2}) = \prod_{i=1}^{N_{\text{bins}}} \text{Poisson}(n_i^{\text{obs}} \mid s_i(f_{M2}) + b_i)
  $$

- Scans $f_{M2}/\Lambda^4$ and finds 95% CL exclusion:

  $$
  q(f_{M2}) = -2 \log \frac{\mathcal{L}(f_{M2})}{\mathcal{L}(0)} > 3.84
  $$

---

## 📊 Output

- `profile_likelihood_fm2_plot.pdf` showing $q(f_{M2})$
- 95% CL bounds printed in console:

```
📉 95% CL Excluded FM2 Range: [-0.3733, +0.3901] TeV⁻⁴
```

---

## 🧾 Citation

If you use this repository, please cite:

> Hamzeh Khanpour, LHeC Semi-leptonic WW Study — Sensitivity to aQGC via BDT + Likelihood Analysis, 2025

---

## 🧪 Requirements

- Python 3.10+
- ROOT with Python bindings
- `xgboost`, `matplotlib`, `scipy`, `pandas`, `numpy`

---

Made with ❤️ by Hamzeh Khanpour
