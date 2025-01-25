#!/usr/bin/env python

import ROOT
from ROOT import TLorentzVector
import numpy as np  # Import numpy for calculations
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


#=========================================================================
#=========================================================================


def process_file(
    file_path,
    hist_lepton,
    hist_leading_jet,
    hist_lepton_eta,
    hist_delta_r=None,
    hist_missing_et=None,
    hist_centrality=None,
    hist_exp_centrality=None,
    hist_jet_centrality=None,
    hist_delta_eta_jj=None,
):
    # Open the ROOT file
    chain = ROOT.TChain("Delphes")
    chain.Add(file_path)

    # Create ExRootTreeReader object
    treeReader = ROOT.ExRootTreeReader(chain)
    numberOfEntries = treeReader.GetEntries()

    # Counters for efficiency calculation
    total_events = numberOfEntries
    selected_events = 0

    # Get branches for electrons, muons, jets, and MET
    branchElectron = treeReader.UseBranch("Electron")
    branchMuon = treeReader.UseBranch("Muon")
    branchJet = treeReader.UseBranch("Jet")
    branchMissingET = treeReader.UseBranch("MissingET")

    # Process each event
    for entry in range(numberOfEntries):
        treeReader.ReadEntry(entry)

        # Count the number of leptons (electrons + muons)
        leptons = []
        for i in range(branchElectron.GetEntries()):
            electron = branchElectron.At(i)
            lepton_vec = TLorentzVector()
            lepton_vec.SetPtEtaPhiM(electron.PT, electron.Eta, electron.Phi, 0.0)
            if electron.PT > 10:
                leptons.append(lepton_vec)
        for i in range(branchMuon.GetEntries()):
            muon = branchMuon.At(i)
            lepton_vec = TLorentzVector()
            lepton_vec.SetPtEtaPhiM(muon.PT, muon.Eta, muon.Phi, 0.0)
            if muon.PT > 10:
                leptons.append(lepton_vec)

        # Count the number of jets
        jets = []
        for i in range(branchJet.GetEntries()):
            jet = branchJet.At(i)
            jet_vec = TLorentzVector()
            jet_vec.SetPtEtaPhiM(jet.PT, jet.Eta, jet.Phi, jet.Mass)
            if jet.PT > 20:
                jets.append(jet_vec)

        # Apply selection criteria: exactly one lepton and exactly two jets
        if len(leptons) != 1 or len(jets) != 2:
            continue

        # Count selected events
        selected_events += 1

        # Fill histograms
        hist_lepton.Fill(leptons[0].Pt())
        hist_lepton_eta.Fill(leptons[0].Eta())
        leading_jet = jets[0] if jets[0].Pt() > jets[1].Pt() else jets[1]
        hist_leading_jet.Fill(leading_jet.Pt())

        # Fill optional histograms
        if hist_delta_r is not None:
            delta_r = leptons[0].DeltaR(leading_jet)
            hist_delta_r.Fill(delta_r)

        if hist_missing_et is not None and branchMissingET.GetEntries() > 0:
            missing_et = branchMissingET.At(0)
            hist_missing_et.Fill(missing_et.MET)

        if hist_centrality is not None and hist_exp_centrality is not None:
            delta_eta_jj = abs(jets[0].Eta() - jets[1].Eta())
            if delta_eta_jj > 0:
                centrality = (leptons[0].Eta() - (jets[0].Eta() + jets[1].Eta()) / 2.0) / delta_eta_jj
                hist_centrality.Fill(centrality)
                hist_exp_centrality.Fill(np.exp(-abs(centrality)))

        if hist_jet_centrality is not None:
            jet_centrality = abs(jets[0].Eta() + jets[1].Eta()) / 2.0
            hist_jet_centrality.Fill(jet_centrality)

        if hist_delta_eta_jj is not None:
            delta_eta_jj = abs(jets[0].Eta() - jets[1].Eta())
            hist_delta_eta_jj.Fill(delta_eta_jj)

    # Calculate selection efficiency
    efficiency = selected_events / total_events if total_events > 0 else 0

    # Return histograms and efficiency
    return hist_lepton, hist_leading_jet, hist_lepton_eta, efficiency




#=========================================================================
#=========================================================================



# Parameters for differential cross-section
signal_cross_section_0 = 1.994  # pb
background_cross_section = 0.0134802  # pb
num_bins = 50
pt_range_lepton = (0, 500)     # Range for lepton pT
pt_range_jet = (0, 500)        # Range for leading jet pT (adjusted for higher jet momenta)
eta_range = (-10, 10)          # Range for pseudorapidity



# Calculate bin width
bin_width_pt_lepton = (pt_range_lepton[1] - pt_range_lepton[0]) / num_bins
bin_width_pt_jet = (pt_range_jet[1] - pt_range_jet[0]) / num_bins
bin_width_eta = (eta_range[1] - eta_range[0]) / num_bins


# Function to calculate differential cross-section
def calculate_dsigma(histogram, total_cross_section, bin_width):
    counts = [histogram.GetBinContent(i) for i in range(1, histogram.GetNbinsX() + 1)]
    bin_edges = [histogram.GetBinLowEdge(i) for i in range(1, histogram.GetNbinsX() + 2)]
    dsigma = [count * (total_cross_section / sum(counts)) / bin_width for count in counts]
    return bin_edges[:-1], dsigma



#=========================================================================
#=========================================================================



# Create histograms for lepton and leading jet pT
hist_lepton_signal_0 = ROOT.TH1F("hist_lepton_signal_0", "Lepton pT Distribution; p_{T} [GeV]; Entries", num_bins, *pt_range_lepton)
hist_lepton_background = ROOT.TH1F("hist_lepton_background", "Lepton pT Distribution; p_{T} [GeV]; Entries", num_bins, *pt_range_lepton)

hist_leading_jet_signal_0 = ROOT.TH1F("hist_leading_jet_signal_0", "Leading Jet pT Distribution; p_{T} [GeV]; Entries", num_bins, *pt_range_jet)
hist_leading_jet_background = ROOT.TH1F("hist_leading_jet_background", "Leading Jet pT Distribution; p_{T} [GeV]; Entries", num_bins, *pt_range_jet)

hist_lepton_eta_signal_0 = ROOT.TH1F("hist_lepton_eta_signal_0", "Lepton Eta Distribution; #eta; Entries", num_bins, *eta_range)
hist_lepton_eta_background = ROOT.TH1F("hist_lepton_eta_background", "Lepton Eta Distribution; #eta; Entries", num_bins, *eta_range)




# Process signal and background files and populate histograms
hist_lepton_signal_0, hist_leading_jet_signal_0, hist_lepton_eta_signal_0, signal_efficiency = process_file(
    signal_file_path, hist_lepton_signal_0, hist_leading_jet_signal_0, hist_lepton_eta_signal_0
)

hist_lepton_background, hist_leading_jet_background, hist_lepton_eta_background, background_efficiency = process_file(
    background_file_path, hist_lepton_background, hist_leading_jet_background, hist_lepton_eta_background
)


# Print selection efficiencies
print(f"Signal Selection Efficiency: {signal_efficiency:.2%}")
print(f"Background Selection Efficiency: {background_efficiency:.2%}")




# Calculate differential cross-sections for lepton pT
pt_bins_lepton_signal_0, dsigma_lepton_signal_0 = calculate_dsigma(
    hist_lepton_signal_0, signal_cross_section_0, bin_width_pt_lepton
)
# Calculate differential cross-sections for leading jet pT
pt_bins_jet_signal_0, dsigma_jet_signal_0 = calculate_dsigma(
    hist_leading_jet_signal_0, signal_cross_section_0, bin_width_pt_jet
)

# Calculate differential cross-sections for lepton η
eta_bins_lepton_signal_0, dsigma_eta_signal_0 = calculate_dsigma(
    hist_lepton_eta_signal_0, signal_cross_section_0, bin_width_eta
)




# Background distributions
pt_bins_lepton_background, dsigma_lepton_background = calculate_dsigma(
    hist_lepton_background, background_cross_section, bin_width_pt_lepton
)
pt_bins_jet_background, dsigma_jet_background = calculate_dsigma(
    hist_leading_jet_background, background_cross_section, bin_width_pt_jet
)

# Calculate differential cross-sections for lepton η (background)
eta_bins_lepton_background, dsigma_eta_background = calculate_dsigma(
    hist_lepton_eta_background, background_cross_section, bin_width_eta
)




#=========================================================================
#=========================================================================



#fig, ax = plt.subplots(figsize=(12.0, 10.0))
plt.subplots_adjust(left=0.15, right=0.95, bottom=0.12, top=0.95)




# Plot the differential cross-sections for lepton pT
#plt.figure(figsize=(8, 9))

plt.step(pt_bins_lepton_signal_0, dsigma_lepton_signal_0, where="mid", alpha=0.7, label="LHeC@1.2 TeV : Signal ($w^+ w^-) [f_{M_0} / \Lambda^4$]", color="red", linewidth=3)
plt.step(pt_bins_lepton_background, dsigma_lepton_background, where="mid", alpha=0.7, label="LHeC@1.2 TeV : SM background ($w^+ w^-$)", color="blue", linewidth=3)
plt.xlabel(r"$p_T^{\ell} \ \mathrm{[GeV]}$")
plt.ylabel(r"$\frac{d\sigma}{dp_T^{\ell}} \ \mathrm{[pb/GeV]}$")
plt.title(r"$e^- p \to e^- w^+ w^- p \to e^- j j \ell \nu_{\ell} p$ : LHeC@1.2 TeV", fontsize=20)
plt.yscale("log")
plt.legend()
plt.grid(True, linestyle="--", alpha=0.6)
plt.tight_layout()
plt.savefig("/home/hamzeh-khanpour/Documents/GitHub/LHeC_Fast_Simulation/Pythia8_Delphes/differential_cross_section_lepton_pt.png", dpi=600)
plt.show()




# Plot the differential cross-sections for leading jet pT
#plt.figure(figsize=(8, 9))

plt.step(pt_bins_jet_signal_0, dsigma_jet_signal_0, where="mid", alpha=0.7, label="LHeC@1.2 TeV : Signal ($w^+ w^-) [f_{M_0} / \Lambda^4$]", color="red", linewidth=3)
plt.step(pt_bins_jet_background, dsigma_jet_background, where="mid", alpha=0.7, label="LHeC@1.2 TeV : SM background ($w^+ w^-$)", color="blue", linewidth=3)
plt.xlabel(r"$p_T^{\mathrm{leading~jet}} \ \mathrm{[GeV]}$")
plt.ylabel(r"$\frac{d\sigma}{dp_T^{\mathrm{leading~jet}}} \ \mathrm{[pb/GeV]}$")
plt.title(r"$e^- p \to e^- w^+ w^- p \to e^- j j \ell \nu_{\ell} p$ : LHeC@1.2 TeV", fontsize=20)
plt.yscale("log")
plt.legend()
plt.grid(True, linestyle="--", alpha=0.6)
plt.tight_layout()
plt.savefig("/home/hamzeh-khanpour/Documents/GitHub/LHeC_Fast_Simulation/Pythia8_Delphes/differential_cross_section_jet_pt.png", dpi=600)
plt.show()





# Plot the differential cross-sections for lepton η
#plt.figure(figsize=(8, 9))  # Create a new figure for the lepton η plot

plt.step(eta_bins_lepton_signal_0, dsigma_eta_signal_0, where="mid", alpha=0.7, label="LHeC@1.2 TeV : Signal ($w^+ w^-) [f_{M_0} / \Lambda^4$)", color="red", linewidth=3)
plt.step(eta_bins_lepton_background, dsigma_eta_background, where="mid", alpha=0.7, label="LHeC@1.2 TeV : SM background ($w^+ w^-$)", color="blue", linewidth=3)
plt.xlabel(r"$\eta^{\ell}$")
plt.ylabel(r"$\frac{d\sigma}{d\eta^{\ell}} \ \mathrm{[pb]}$")
plt.title(r"$e^- p \to e^- w^+ w^- p \to e^- j j \ell \nu_{\ell} p$ : LHeC@1.2 TeV",fontsize=20)
plt.yscale("log")
plt.legend()
plt.grid(True, linestyle="--", alpha=0.6)
plt.tight_layout()
plt.savefig("/home/hamzeh-khanpour/Documents/GitHub/LHeC_Fast_Simulation/Pythia8_Delphes/differential_cross_section_lepton_eta.png", dpi=600)
plt.show()



#=========================================================================
#=========================================================================



