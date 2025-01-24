#!/usr/bin/env python

import ROOT
from ROOT import TLorentzVector
import matplotlib.pyplot as plt



# Matplotlib configuration for publication-quality plots
import mplhep as hep

hep.style.use("CMS")
#plt.style.use(hep.style.ROOT)

#plt.rcParams["axes.linewidth"] = 1.8
#plt.rcParams["xtick.major.width"] = 1.8
#plt.rcParams["xtick.minor.width"] = 1.8
#plt.rcParams["ytick.major.width"] = 1.8
#plt.rcParams["ytick.minor.width"] = 1.8

#plt.rcParams["xtick.direction"] = "in"
#plt.rcParams["ytick.direction"] = "in"

#plt.rcParams["xtick.labelsize"] = 15
#plt.rcParams["ytick.labelsize"] = 15

#plt.rcParams["legend.fontsize"] = 15

#plt.rcParams['legend.title_fontsize'] = 'x-large'





# Path to the ROOT files
signal_file_path = "aa_ww_semi_leptonic_NP_FM0.root"
background_file_path = "aa_ww_semi_leptonic_SM.root"

# Load Delphes library
ROOT.gSystem.Load("libDelphes")
ROOT.gInterpreter.Declare('#include "classes/DelphesClasses.h"')
ROOT.gInterpreter.Declare('#include "external/ExRootAnalysis/ExRootTreeReader.h"')

# Function to process a ROOT file and return a filled histogram
def process_file(file_path, hist):
    # Open the ROOT file
    chain = ROOT.TChain("Delphes")
    chain.Add(file_path)

    # Create ExRootTreeReader object
    treeReader = ROOT.ExRootTreeReader(chain)
    numberOfEntries = treeReader.GetEntries()

    # Get branches for electrons and muons
    branchElectron = treeReader.UseBranch("Electron")
    branchMuon = treeReader.UseBranch("Muon")

    # Process each event
    for entry in range(numberOfEntries):
        treeReader.ReadEntry(entry)

        # Fill histogram with electron pT
        for i in range(branchElectron.GetEntries()):
            electron = branchElectron.At(i)
            electron_vec = TLorentzVector()
            electron_vec.SetPtEtaPhiM(electron.PT, electron.Eta, electron.Phi, 0.0)
            hist.Fill(electron_vec.Pt())

        # Fill histogram with muon pT
        for i in range(branchMuon.GetEntries()):
            muon = branchMuon.At(i)
            muon_vec = TLorentzVector()
            muon_vec.SetPtEtaPhiM(muon.PT, muon.Eta, muon.Phi, 0.0)
            hist.Fill(muon_vec.Pt())

    return hist

# Create histograms for signal and background
hist_signal = ROOT.TH1F("hist_signal", "Lepton pT Distribution; p_{T} [GeV]; Entries", 50, 0, 200)
hist_background = ROOT.TH1F("hist_background", "Lepton pT Distribution; p_{T} [GeV]; Entries", 50, 0, 200)

# Process signal and background files
process_file(signal_file_path, hist_signal)
process_file(background_file_path, hist_background)

# Convert ROOT histograms to matplotlib
x_vals_signal = [hist_signal.GetBinCenter(i) for i in range(1, hist_signal.GetNbinsX() + 1)]
y_vals_signal = [hist_signal.GetBinContent(i) for i in range(1, hist_signal.GetNbinsX() + 1)]

x_vals_background = [hist_background.GetBinCenter(i) for i in range(1, hist_background.GetNbinsX() + 1)]
y_vals_background = [hist_background.GetBinContent(i) for i in range(1, hist_background.GetNbinsX() + 1)]

# Plot the histograms
plt.subplots_adjust(left=0.15, right=0.95, bottom=0.12, top=0.95)
plt.figure(figsize=(8, 9))  # Adjust plot size to 10x12

plt.bar(
    x_vals_signal,
    y_vals_signal,
    width=hist_signal.GetBinWidth(1),
    alpha=0.6,
    label="LHeC@1.2 TeV : Signal ($w^+ w^-) [f_{M_0} / \Lambda^4$]", color="red", linewidth=3,
)
plt.bar(
    x_vals_background,
    y_vals_background,
    width=hist_background.GetBinWidth(1),
    alpha=0.6,
    label="LHeC@1.2 TeV : SM background ($w^+ w^-$)", color="blue", linewidth=3,
)

# Add labels, legend, and title
plt.xlabel(r"$p_T^{\ell} \ \mathrm{[GeV]}$")
plt.ylabel("Entries")
plt.title(r"$e^- p \to e^- w^+ w^- p \to e^- j j \ell \nu_{\ell} p$ : LHeC@1.2 TeV", fontsize=24)
plt.legend()
plt.grid()

# Save the plot as a PNG file
output_file = "/home/hamzeh-khanpour/Documents/GitHub/LHeC_Fast_Simulation/Pythia8_Delphes/lepton_pt_comparison.png"
plt.savefig(output_file, dpi=300)
print(f"Plot saved as {output_file}")

# Show the plot (optional)
plt.show()
