#!/usr/bin/env python

import ROOT
from ROOT import TLorentzVector
import matplotlib.pyplot as plt

# Matplotlib configuration for publication-quality plots
import mplhep as hep

hep.style.use("CMS")

# Path to the ROOT files
signal_file_path = "aa_ww_semi_leptonic_NP_FM0.root"
background_file_path = "aa_ww_semi_leptonic_SM.root"





# Load Delphes library
ROOT.gSystem.Load("libDelphes")
ROOT.gInterpreter.Declare('#include "classes/DelphesClasses.h"')
ROOT.gInterpreter.Declare('#include "external/ExRootAnalysis/ExRootTreeReader.h"')



# Function to process a ROOT file and fill histograms for leptons and leading jet
def process_file(file_path, hist_lepton, hist_leading_jet):
    # Open the ROOT file
    chain = ROOT.TChain("Delphes")
    chain.Add(file_path)

    # Create ExRootTreeReader object
    treeReader = ROOT.ExRootTreeReader(chain)
    numberOfEntries = treeReader.GetEntries()

    # Counters for efficiency calculation
    total_events = numberOfEntries
    selected_events = 0

    # Get branches for electrons, muons, and jets
    branchElectron = treeReader.UseBranch("Electron")
    branchMuon = treeReader.UseBranch("Muon")
    branchJet = treeReader.UseBranch("Jet")

    # Process each event
    for entry in range(numberOfEntries):
        treeReader.ReadEntry(entry)

        # Count the number of leptons (electrons + muons)
        leptons = []
        for i in range(branchElectron.GetEntries()):
            electron = branchElectron.At(i)
            lepton_vec = TLorentzVector()
            lepton_vec.SetPtEtaPhiM(electron.PT, electron.Eta, electron.Phi, 0.0)
            if electron.PT > 10:  # Corrected syntax
                leptons.append(lepton_vec)
        for i in range(branchMuon.GetEntries()):
            muon = branchMuon.At(i)
            lepton_vec = TLorentzVector()
            lepton_vec.SetPtEtaPhiM(muon.PT, muon.Eta, muon.Phi, 0.0)
            if muon.PT > 10:  # Corrected syntax
                leptons.append(lepton_vec)

        # Count the number of jets
        jets = []
        for i in range(branchJet.GetEntries()):
            jet = branchJet.At(i)
            jet_vec = TLorentzVector()
            jet_vec.SetPtEtaPhiM(jet.PT, jet.Eta, jet.Phi, jet.Mass)
            if jet.PT > 20:  # Added jet pT threshold for selection
                jets.append(jet_vec)

        # Apply selection criteria: exactly one lepton and exactly two jets
        if len(leptons) != 1 or len(jets) != 2:
            continue  # Skip this event if it doesn't meet the criteria

        # Count selected events
        selected_events += 1

        # Fill histogram with lepton pT
        hist_lepton.Fill(leptons[0].Pt())

        # Fill histogram with leading jet pT
        leading_jet = jets[0] if jets[0].Pt() > jets[1].Pt() else jets[1]
        hist_leading_jet.Fill(leading_jet.Pt())

    # Calculate selection efficiency
    efficiency = selected_events / total_events if total_events > 0 else 0

    return hist_lepton, hist_leading_jet, efficiency







# Create histograms for lepton and leading jet pT (signal and background)
hist_lepton_signal = ROOT.TH1F("hist_lepton_signal", "Lepton pT Distribution; p_{T} [GeV]; Entries", 50, 0, 400)
hist_lepton_background = ROOT.TH1F("hist_lepton_background", "Lepton pT Distribution; p_{T} [GeV]; Entries", 50, 0, 400)


hist_leading_jet_signal = ROOT.TH1F("hist_leading_jet_signal", "Leading Jet pT Distribution; p_{T} [GeV]; Entries", 50, 0, 400)
hist_leading_jet_background = ROOT.TH1F("hist_leading_jet_background", "Leading Jet pT Distribution; p_{T} [GeV]; Entries", 50, 0, 400)


# Process signal and background files and calculate efficiencies
hist_lepton_signal, hist_leading_jet_signal, signal_efficiency = process_file(signal_file_path, hist_lepton_signal, hist_leading_jet_signal)
hist_lepton_background, hist_leading_jet_background, background_efficiency = process_file(background_file_path, hist_lepton_background, hist_leading_jet_background)



# Print selection efficiencies
print(f"Signal Selection Efficiency: {signal_efficiency:.2%}")
print(f"Background Selection Efficiency: {background_efficiency:.2%}")



# Convert ROOT histograms to matplotlib for leptons
x_vals_lepton_signal = [hist_lepton_signal.GetBinCenter(i) for i in range(1, hist_lepton_signal.GetNbinsX() + 1)]
y_vals_lepton_signal = [hist_lepton_signal.GetBinContent(i) for i in range(1, hist_lepton_signal.GetNbinsX() + 1)]

x_vals_lepton_background = [hist_lepton_background.GetBinCenter(i) for i in range(1, hist_lepton_background.GetNbinsX() + 1)]
y_vals_lepton_background = [hist_lepton_background.GetBinContent(i) for i in range(1, hist_lepton_background.GetNbinsX() + 1)]

# Convert ROOT histograms to matplotlib for leading jet
x_vals_jet_signal = [hist_leading_jet_signal.GetBinCenter(i) for i in range(1, hist_leading_jet_signal.GetNbinsX() + 1)]
y_vals_jet_signal = [hist_leading_jet_signal.GetBinContent(i) for i in range(1, hist_leading_jet_signal.GetNbinsX() + 1)]

x_vals_jet_background = [hist_leading_jet_background.GetBinCenter(i) for i in range(1, hist_leading_jet_background.GetNbinsX() + 1)]
y_vals_jet_background = [hist_leading_jet_background.GetBinContent(i) for i in range(1, hist_leading_jet_background.GetNbinsX() + 1)]





# Plot the histograms for leptons
plt.figure(figsize=(8, 9))  # Create a new figure for the lepton plot
plt.subplots_adjust(left=0.15, right=0.95, bottom=0.12, top=0.95)

plt.bar(
    x_vals_lepton_signal,
    y_vals_lepton_signal,
    width=hist_lepton_signal.GetBinWidth(1),
    alpha=0.6,
    label="LHeC@1.2 TeV : Signal ($w^+ w^-) [f_{M_0} / \Lambda^4$)", color="red",
)
plt.bar(
    x_vals_lepton_background,
    y_vals_lepton_background,
    width=hist_lepton_background.GetBinWidth(1),
    alpha=0.6,
    label="LHeC@1.2 TeV : SM background ($w^+ w^-$)", color="blue",
)

plt.xlabel(r"$p_T^{\ell} \ \mathrm{[GeV]}$")
plt.ylabel("Entries")
plt.title(r"$e^- p \to e^- w^+ w^- p \to e^- j j \ell \nu_{\ell} p$ : LHeC@1.2 TeV", fontsize=18)
plt.legend()
plt.grid()
plt.savefig("/home/hamzeh-khanpour/Documents/GitHub/LHeC_Fast_Simulation/Pythia8_Delphes/lepton_pt_comparison.png", dpi=300)
plt.show()  # Ensure this ends the current plot properly




# Plot the histograms for leading jets
plt.figure(figsize=(8, 9))  # Create a new figure for the jet plot
plt.subplots_adjust(left=0.15, right=0.95, bottom=0.12, top=0.95)

plt.bar(
    x_vals_jet_signal,
    y_vals_jet_signal,
    width=hist_leading_jet_signal.GetBinWidth(1),
    alpha=0.6,
    label="LHeC@1.2 TeV : Signal ($w^+ w^-) [f_{M_0} / \Lambda^4$)", color="red",
)
plt.bar(
    x_vals_jet_background,
    y_vals_jet_background,
    width=hist_leading_jet_background.GetBinWidth(1),
    alpha=0.6,
    label="LHeC@1.2 TeV : SM background ($w^+ w^-$)", color="blue",
)

plt.xlabel(r"$p_T^{\mathrm{leading~jet}} \ \mathrm{[GeV]}$")
plt.ylabel("Entries")
plt.title(r"$e^- p \to e^- w^+ w^- p \to e^- j j \ell \nu_{\ell} p$ : LHeC@1.2 TeV", fontsize=18)
plt.legend()
plt.grid()
plt.savefig("/home/hamzeh-khanpour/Documents/GitHub/LHeC_Fast_Simulation/Pythia8_Delphes/leading_jet_pt_comparison.png", dpi=300)
plt.show()  # Ensure the second plot is displayed


