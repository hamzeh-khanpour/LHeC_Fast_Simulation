#!/usr/bin/env python

import ROOT
from ROOT import TLorentzVector
import matplotlib.pyplot as plt


# Path to the ROOT file
file_path = "aa_ww_semi_leptonic_SM.root"


# Load Delphes library
ROOT.gSystem.Load("libDelphes")
ROOT.gInterpreter.Declare('#include "classes/DelphesClasses.h"')
ROOT.gInterpreter.Declare('#include "external/ExRootAnalysis/ExRootTreeReader.h"')


# Open the ROOT file and get the Delphes tree
root_file = ROOT.TFile.Open(file_path)
if not root_file or root_file.IsZombie():
    print(f"Error: Cannot open file {file_path}")
    exit(1)

chain = ROOT.TChain("Delphes")
chain.Add(file_path)


# Create an ExRootTreeReader object
treeReader = ROOT.ExRootTreeReader(chain)
numberOfEntries = treeReader.GetEntries()



# Get branches for electrons and muons
branchElectron = treeReader.UseBranch("Electron")
branchMuon = treeReader.UseBranch("Muon")


# Create a histogram for lepton transverse momentum (pT)
h_pt = ROOT.TH1F("h_pt", "Lepton Transverse Momentum; p_{T} [GeV]; Entries", 50, 0, 200)


# Loop over all events
for entry in range(numberOfEntries):
    # Load the data for the current event
    treeReader.ReadEntry(entry)

    # Loop over all electrons in the event
    for i in range(branchElectron.GetEntries()):
        electron = branchElectron.At(i)

        # Create a TLorentzVector for the electron
        electron_vec = TLorentzVector()
        electron_vec.SetPtEtaPhiM(electron.PT, electron.Eta, electron.Phi, 0.0)

        # Fill the histogram with pT
        h_pt.Fill(electron_vec.Pt())

        # Print electron properties (optional for debugging)
        print(f"Electron: pT = {electron_vec.Pt():.2f}, Eta = {electron_vec.Eta():.2f}, Phi = {electron_vec.Phi():.2f}, Mass = {electron_vec.M():.2f}")

    # Loop over all muons in the event
    for i in range(branchMuon.GetEntries()):
        muon = branchMuon.At(i)

        # Create a TLorentzVector for the muon
        muon_vec = TLorentzVector()
        muon_vec.SetPtEtaPhiM(muon.PT, muon.Eta, muon.Phi, 0.0)

        # Fill the histogram with pT
        h_pt.Fill(muon_vec.Pt())

        # Print muon properties (optional for debugging)
#        print(f"Muon: pT = {muon_vec.Pt():.2f}, Eta = {muon_vec.Eta():.2f}, Phi = {muon_vec.Phi():.2f}, Mass = {muon_vec.M():.2f}")



# Convert the ROOT histogram to matplotlib for visualization
x_vals = [h_pt.GetBinCenter(i) for i in range(1, h_pt.GetNbinsX() + 1)]
y_vals = [h_pt.GetBinContent(i) for i in range(1, h_pt.GetNbinsX() + 1)]


# Plot the histogram
plt.bar(x_vals, y_vals, width=h_pt.GetBinWidth(1), color='blue', alpha=0.7, label="Leptons")
plt.xlabel("Transverse Momentum (p_{T}) [GeV]")
plt.ylabel("Entries")
plt.title("Lepton pT Distribution")
plt.legend()
plt.grid()
plt.show()

# Close the ROOT file
root_file.Close()

