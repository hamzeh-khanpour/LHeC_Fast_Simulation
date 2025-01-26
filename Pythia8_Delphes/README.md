# LHeC Fast Simulation: Semi-Leptonic Decay of \( W^+ W^- \)

This repository provides an analysis framework for simulating and analyzing semi-leptonic decays of \( W^+ W^- \) at the Large Hadron Electron Collider (LHeC). The project uses **Delphes** for fast detector simulation and **ROOT** for data analysis and histogramming.

---

## Table of Contents

- [Overview](#overview)
- [Features](#features)
- [Prerequisites](#prerequisites)
- [Setup](#setup)
- [Usage](#usage)
- [Key Distributions](#key-distributions)
- [Output](#output)
- [Directory Structure](#directory-structure)
- [Acknowledgments](#acknowledgments)

---

## Overview

This project focuses on the semi-leptonic decay channel of \( W^+ W^- \), where:
- One \( W \)-boson decays hadronically (\( W \to jj \)).
- The other \( W \)-boson decays leptonically (\( W \to \ell \nu \)).

**Goals**:
1. Perform event selection based on specific criteria.
2. Compute differential cross-sections for various observables.
3. Save results as ROOT histograms and publication-ready plots.

---

## Features

### Selection Criteria
- Events must satisfy the following:
  - **1 lepton** (electron or muon) with \( p_T > 10 \, \text{GeV} \).
  - **2 jets** with \( p_T > 20 \, \text{GeV} \).

### Distributions
The framework calculates and plots the following:
1. **Lepton Observables**:
   - Transverse momentum (\( p_T^{\ell} \)).
   - Pseudorapidity (\( \eta^{\ell} \)).

2. **Jet Observables**:
   - Leading jet transverse momentum (\( p_T^{\mathrm{leading~jet}} \)).
   - Jet centrality (\( C_{\mathrm{jets}} \)).

3. **Other Observables**:
   - Missing Transverse Energy (MET).
   - Centrality (\( C_{\ell} \)).
   - Exponential Centrality (\( C_{\ell}^{\mathrm{exp}} \)).
   - \( \Delta R \) between lepton and leading jet.
   - \( \Delta \eta_{jj} \): Pseudorapidity difference between jets.


## Mathematical Formulas

### 1. Lepton Centrality (\( C_{\ell} \))
The lepton centrality measures the pseudorapidity position of the lepton relative to the two jets:
\[
C_{\ell} = \frac{\eta_{\ell} - \frac{\eta_{\mathrm{jet1}} + \eta_{\mathrm{jet2}}}{2}}{\Delta \eta_{jj}}
\]
Where:
- \( \eta_{\ell} \): Pseudorapidity of the lepton.
- \( \eta_{\mathrm{jet1}} \): Pseudorapidity of the first jet.
- \( \eta_{\mathrm{jet2}} \): Pseudorapidity of the second jet.
- \( \Delta \eta_{jj} \): Pseudorapidity difference between the two jets:
  \[
  \Delta \eta_{jj} = |\eta_{\mathrm{jet1}} - \eta_{\mathrm{jet2}}|
  \]

---

### 2. Exponential Centrality (\( C_{\ell}^{\mathrm{exp}} \))
The exponential centrality highlights events with small \( C_{\ell} \) values:
\[
C_{\ell}^{\mathrm{exp}} = e^{-|C_{\ell}|}
\]
Where:
- \( |C_{\ell}| \): Absolute value of the lepton centrality.

---

### 3. Jet Centrality (\( C_{\mathrm{jets}} \))
The jet centrality quantifies the average pseudorapidity position of the two jets:
\[
C_{\mathrm{jets}} = \frac{|\eta_{\mathrm{jet1}} + \eta_{\mathrm{jet2}}|}{2}
\]
Where:
- \( \eta_{\mathrm{jet1}} \): Pseudorapidity of the first jet.
- \( \eta_{\mathrm{jet2}} \): Pseudorapidity of the second jet.




### Histogram Export
- All histograms are saved in `output_histograms.root` with branches for signal and background.

### Mathematical Formula Embedding
- Observables' formulas are included directly in the plots.

---

## Prerequisites

### Software
- **Python 3.10+**
- **ROOT** (with Python bindings)
- **Delphes 3.5.0**


### Python Libraries
Install required packages:
```bash
pip install numpy matplotlib mplhep
