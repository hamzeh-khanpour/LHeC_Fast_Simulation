# LHeC Fast Simulation: Semi-Leptonic Decay of \( W^+ W^- \)

This repository contains a Python-based analysis framework for simulating and analyzing semi-leptonic decay events of \( W^+ W^- \) at the Large Hadron Electron Collider (LHeC). The simulation uses **Delphes** for fast detector simulation and **ROOT** for data processing and histogramming.

---

## Table of Contents

- [Overview](#overview)
- [Features](#features)
- [Prerequisites](#prerequisites)
- [Setup](#setup)
- [Usage](#usage)
- [Key Distributions](#key-distributions)
- [Output](#output)
- [Acknowledgments](#acknowledgments)

---

## Overview

This project analyzes semi-leptonic decays of \( W^+ W^- \), where:
- One \( W \)-boson decays into jets (hadronic decay).
- The other \( W \)-boson decays into a lepton and a neutrino (leptonic decay).

The repository includes functionality to:
- Parse ROOT files produced by Delphes.
- Apply selection criteria to events.
- Calculate differential cross-sections.
- Plot and save histograms for key physical observables.

---

## Features

1. **Selection Criteria**:
   - Events are selected if they have exactly:
     - **1 lepton** (electron or muon) with \( p_T > 10 \, \text{GeV} \).
     - **2 jets** with \( p_T > 20 \, \text{GeV} \).

2. **Distributions**:
   - Lepton \( p_T \) and \( \eta \).
   - Leading jet \( p_T \).
   - Missing Transverse Energy (MET).
   - Centrality (\( C_{\ell} \)).
   - Exponential Centrality (\( C_{\ell}^{\mathrm{exp}} \)).
   - Jet Centrality (\( C_{\mathrm{jets}} \)).
   - \( \Delta R \) (distance between lepton and leading jet).
   - \( \Delta \eta_{jj} \) (pseudorapidity difference between jets).

3. **ROOT Histogram Export**:
   - All histograms are saved into `output_histograms.root` with branches for signal and background.

4. **Mathematical Formulas in Plots**:
   - Mathematical definitions of observables are embedded directly in the plots.

---

## Prerequisites

- **Python 3.10+**
- **ROOT** (with Python bindings)
- **Delphes 3.5.0**
- Python Libraries:
  - `numpy`
  - `matplotlib`
  - `mplhep`

---

## Setup

1. Clone the repository:
   ```bash
   git clone https://github.com/your_username/LHeC_Fast_Simulation.git
   cd LHeC_Fast_Simulation
