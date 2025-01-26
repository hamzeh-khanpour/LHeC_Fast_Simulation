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
signal_file_path = "aa_ww_semi_leptonic_NP_FM0_new.root"
background_file_path = "aa_ww_semi_leptonic_SM_new.root"



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
    hist_delta_r,
    hist_missing_et,
    hist_centrality,
    hist_exp_centrality,
    hist_jet_centrality,
    hist_delta_eta_jj,
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
            if jet.PT > 10:
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
    return hist_lepton, hist_leading_jet, hist_lepton_eta, hist_delta_r, hist_missing_et, hist_centrality, hist_exp_centrality, hist_jet_centrality, hist_delta_eta_jj, efficiency





#=========================================================================
#=========================================================================



# Parameters for differential cross-section
signal_cross_section_0   = 1.994      # pb
background_cross_section = 0.0134802  # pb
num_bins = 50
pt_range_lepton = (0, 500)     # Range for lepton pT
pt_range_jet = (0, 500)        # Range for leading jet pT (adjusted for higher jet momenta)
eta_range = (-10, 10)          # Range for pseudorapidity
delta_r_range = (0, 10)        # Range for Delta R
met_range = (0, 500)           # Range for Missing Transverse Energy (MET)
centrality_range = (-5, 5)     # Range for centrality
exp_centrality_range = (0, 2)  # Range for exponential centrality
jet_centrality_range = (0, 5)  # Range for jet centrality
delta_eta_jj_range = (0, 6)    # Range for Delta Eta between jets




# Calculate bin width
bin_width_pt_lepton = (pt_range_lepton[1] - pt_range_lepton[0]) / num_bins
bin_width_pt_jet = (pt_range_jet[1] - pt_range_jet[0]) / num_bins
bin_width_eta = (eta_range[1] - eta_range[0]) / num_bins
bin_width_delta_r = (delta_r_range[1] - delta_r_range[0]) / num_bins
bin_width_met = (met_range[1] - met_range[0]) / num_bins
bin_width_centrality = (centrality_range[1] - centrality_range[0]) / num_bins
bin_width_exp_centrality = (exp_centrality_range[1] - exp_centrality_range[0]) / num_bins
bin_width_jet_centrality = (jet_centrality_range[1] - jet_centrality_range[0]) / num_bins
bin_width_delta_eta_jj = (delta_eta_jj_range[1] - delta_eta_jj_range[0]) / num_bins





# Function to calculate differential cross-section
def calculate_dsigma(histogram, total_cross_section, bin_width):
    counts = [histogram.GetBinContent(i) for i in range(1, histogram.GetNbinsX() + 1)]
    bin_edges = [histogram.GetBinLowEdge(i) for i in range(1, histogram.GetNbinsX() + 2)]

    if sum(counts) == 0:  # Avoid division by zero
        print("Warning: Histogram is empty. Returning zero differential cross-section.")
        dsigma = [0] * len(counts)
    else:
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


# Create histograms for Delta R
hist_delta_r_signal_0 = ROOT.TH1F("hist_delta_r_signal_0", "Delta R Distribution (Signal); ΔR; Entries", num_bins, *delta_r_range)
hist_delta_r_background = ROOT.TH1F("hist_delta_r_background", "Delta R Distribution (Background); ΔR; Entries", num_bins, *delta_r_range)


# Create histograms for Missing Transverse Energy (MET)
hist_missing_et_signal_0 = ROOT.TH1F("hist_missing_et_signal_0", "Missing ET Distribution (Signal); MET [GeV]; Entries", num_bins, *met_range)
hist_missing_et_background = ROOT.TH1F("hist_missing_et_background", "Missing ET Distribution (Background); MET [GeV]; Entries", num_bins, *met_range)

# Create histograms for centrality
hist_centrality_signal_0 = ROOT.TH1F("hist_centrality_signal_0", "Centrality Distribution (Signal); Centrality; Entries", num_bins, *centrality_range)
hist_centrality_background = ROOT.TH1F("hist_centrality_background", "Centrality Distribution (Background); Centrality; Entries", num_bins, *centrality_range)


# Create histograms for exponential centrality
hist_exp_centrality_signal_0 = ROOT.TH1F("hist_exp_centrality_signal_0", "Exponential Centrality (Signal); exp(Centrality); Entries", num_bins, *exp_centrality_range)
hist_exp_centrality_background = ROOT.TH1F("hist_exp_centrality_background", "Exponential Centrality (Background); exp(Centrality); Entries", num_bins, *exp_centrality_range)


# Create histograms for jet centrality
hist_jet_centrality_signal_0 = ROOT.TH1F("hist_jet_centrality_signal_0", "Jet Centrality Distribution (Signal); Jet Centrality; Entries", num_bins, *jet_centrality_range)
hist_jet_centrality_background = ROOT.TH1F("hist_jet_centrality_background", "Jet Centrality Distribution (Background); Jet Centrality; Entries", num_bins, *jet_centrality_range)


# Create histograms for Delta Eta between jets (Δηjj)
hist_delta_eta_jj_signal_0 = ROOT.TH1F("hist_delta_eta_jj_signal_0", "Delta Eta Between Jets (Signal); Δηjj; Entries", num_bins, *delta_eta_jj_range)
hist_delta_eta_jj_background = ROOT.TH1F("hist_delta_eta_jj_background", "Delta Eta Between Jets (Background); Δηjj; Entries", num_bins, *delta_eta_jj_range)





# Process signal and background files and populate histograms
(hist_lepton_signal_0, hist_leading_jet_signal_0, hist_lepton_eta_signal_0, hist_delta_r_signal_0,
 hist_missing_et_signal_0, hist_centrality_signal_0, hist_exp_centrality_signal_0, hist_jet_centrality_signal_0,
 hist_delta_eta_jj_signal_0, signal_efficiency) = process_file(
    signal_file_path, hist_lepton_signal_0, hist_leading_jet_signal_0, hist_lepton_eta_signal_0,
    hist_delta_r_signal_0, hist_missing_et_signal_0, hist_centrality_signal_0, hist_exp_centrality_signal_0,
    hist_jet_centrality_signal_0, hist_delta_eta_jj_signal_0
)


(hist_lepton_background, hist_leading_jet_background, hist_lepton_eta_background, hist_delta_r_background,
 hist_missing_et_background, hist_centrality_background, hist_exp_centrality_background, hist_jet_centrality_background,
 hist_delta_eta_jj_background, background_efficiency) = process_file(
    background_file_path, hist_lepton_background, hist_leading_jet_background, hist_lepton_eta_background,
    hist_delta_r_background, hist_missing_et_background, hist_centrality_background, hist_exp_centrality_background,
    hist_jet_centrality_background, hist_delta_eta_jj_background
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

# Calculate differential cross-sections for Delta R
delta_r_bins_signal_0, dsigma_delta_r_signal_0 = calculate_dsigma(
    hist_delta_r_signal_0, signal_cross_section_0, bin_width_delta_r
)

# Calculate differential cross-sections for MET
met_bins_signal_0, dsigma_met_signal_0 = calculate_dsigma(
    hist_missing_et_signal_0, signal_cross_section_0, bin_width_met
)

# Calculate differential cross-sections for centrality
centrality_bins_signal_0, dsigma_centrality_signal_0 = calculate_dsigma(
    hist_centrality_signal_0, signal_cross_section_0, bin_width_centrality
)

# Calculate differential cross-sections for exponential centrality
exp_centrality_bins_signal_0, dsigma_exp_centrality_signal_0 = calculate_dsigma(
    hist_exp_centrality_signal_0, signal_cross_section_0, bin_width_exp_centrality
)

# Calculate differential cross-sections for jet centrality
jet_centrality_bins_signal_0, dsigma_jet_centrality_signal_0 = calculate_dsigma(
    hist_jet_centrality_signal_0, signal_cross_section_0, bin_width_jet_centrality
)

# Calculate differential cross-sections for Δηjj
delta_eta_jj_bins_signal_0, dsigma_delta_eta_jj_signal_0 = calculate_dsigma(
    hist_delta_eta_jj_signal_0, signal_cross_section_0, bin_width_delta_eta_jj
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

delta_r_bins_background, dsigma_delta_r_background = calculate_dsigma(
    hist_delta_r_background, background_cross_section, bin_width_delta_r
)

met_bins_background, dsigma_met_background = calculate_dsigma(
    hist_missing_et_background, background_cross_section, bin_width_met
)

centrality_bins_background, dsigma_centrality_background = calculate_dsigma(
    hist_centrality_background, background_cross_section, bin_width_centrality
)

exp_centrality_bins_background, dsigma_exp_centrality_background = calculate_dsigma(
    hist_exp_centrality_background, background_cross_section, bin_width_exp_centrality
)

jet_centrality_bins_background, dsigma_jet_centrality_background = calculate_dsigma(
    hist_jet_centrality_background, background_cross_section, bin_width_jet_centrality
)

delta_eta_jj_bins_background, dsigma_delta_eta_jj_background = calculate_dsigma(
    hist_delta_eta_jj_background, background_cross_section, bin_width_delta_eta_jj
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
plt.title(r"Delphes simulation : $e^- p \to e^- w^+ w^- p \to e^- j j \ell \nu_{\ell} p$ : LHeC@1.2 TeV", fontsize=20)
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
plt.title(r"Delphes simulation : $e^- p \to e^- w^+ w^- p \to e^- j j \ell \nu_{\ell} p$ : LHeC@1.2 TeV", fontsize=20)
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
plt.title(r"Delphes simulation : $e^- p \to e^- w^+ w^- p \to e^- j j \ell \nu_{\ell} p$ : LHeC@1.2 TeV",fontsize=20)
plt.yscale("log")
plt.legend()
plt.grid(True, linestyle="--", alpha=0.6)
plt.tight_layout()
plt.savefig("/home/hamzeh-khanpour/Documents/GitHub/LHeC_Fast_Simulation/Pythia8_Delphes/differential_cross_section_lepton_eta.png", dpi=600)
plt.show()





# Plot the differential cross-sections for Delta R
#plt.figure(figsize=(8, 9))

plt.step(delta_r_bins_signal_0, dsigma_delta_r_signal_0, where="mid", alpha=0.7, label="LHeC@1.2 TeV : Signal ($w^+ w^-) [f_{M_0} / \Lambda^4$)", color="red", linewidth=3)
plt.step(delta_r_bins_background, dsigma_delta_r_background, where="mid", alpha=0.7, label="LHeC@1.2 TeV : SM background ($w^+ w^-$)", color="blue", linewidth=3)
plt.xlabel(r"$\Delta R(\ell, \mathrm{leading~jet})$")
plt.ylabel(r"$\frac{d\sigma}{d\Delta R} \ \mathrm{[pb]}$")
plt.title(r"Delphes simulation : $e^- p \to e^- w^+ w^- p \to e^- j j \ell \nu_{\ell} p$ : LHeC@1.2 TeV", fontsize=20)
plt.yscale("log")
plt.legend()
plt.grid(True, linestyle="--", alpha=0.6)
plt.tight_layout()
plt.savefig("/home/hamzeh-khanpour/Documents/GitHub/LHeC_Fast_Simulation/Pythia8_Delphes/differential_cross_section_delta_r.png", dpi=600)
plt.show()




# Plot the differential cross-sections for MET
#plt.figure(figsize=(8, 9))

plt.step(met_bins_signal_0, dsigma_met_signal_0, where="mid", alpha=0.7, label="LHeC@1.2 TeV : Signal ($w^+ w^-) [f_{M_0} / \Lambda^4$)", color="red", linewidth=3)
plt.step(met_bins_background, dsigma_met_background, where="mid", alpha=0.7, label="LHeC@1.2 TeV : SM background ($w^+ w^-$)", color="blue", linewidth=3)
plt.xlabel(r"$\mathrm{MET} \ \mathrm{[GeV]}$")
plt.ylabel(r"$\frac{d\sigma}{d\mathrm{MET}} \ \mathrm{[pb/GeV]}$")
plt.title(r"Delphes simulation : $e^- p \to e^- w^+ w^- p \to e^- j j \ell \nu_{\ell} p$ : LHeC@1.2 TeV", fontsize=20)
plt.yscale("log")
plt.legend()
plt.grid(True, linestyle="--", alpha=0.6)
plt.tight_layout()
plt.savefig("/home/hamzeh-khanpour/Documents/GitHub/LHeC_Fast_Simulation/Pythia8_Delphes/differential_cross_section_met.png", dpi=600)
plt.show()






# Plot the differential cross-sections for centrality
#plt.figure(figsize=(8, 9))

plt.step(centrality_bins_signal_0, dsigma_centrality_signal_0, where="mid", alpha=0.7, label="LHeC@1.2 TeV : Signal ($w^+ w^-) [f_{M_0} / \Lambda^4$)", color="red", linewidth=3)
plt.step(centrality_bins_background, dsigma_centrality_background, where="mid", alpha=0.7, label="LHeC@1.2 TeV : SM background ($w^+ w^-$)", color="blue", linewidth=3)
plt.xlabel(r"$C_{\ell}$")
plt.ylabel(r"$\frac{d\sigma}{dC_{\ell}} \ \mathrm{[pb]}$")
plt.title(r"Delphes simulation : $e^- p \to e^- w^+ w^- p \to e^- j j \ell \nu_{\ell} p$ : LHeC@1.2 TeV", fontsize=20)
plt.yscale("log")
plt.legend()
plt.grid(True, linestyle="--", alpha=0.6)
# Add formula inside the plot
plt.text(0.5, 0.001, r"$C_{\ell} = \frac{\eta_{\ell} - \frac{\eta_{\mathrm{jet1}} + \eta_{\mathrm{jet2}}}{2}}{\Delta \eta_{jj}}$",  color="black")
plt.tight_layout()
plt.savefig("/home/hamzeh-khanpour/Documents/GitHub/LHeC_Fast_Simulation/Pythia8_Delphes/differential_cross_section_centrality.png", dpi=600)
plt.show()





# Plot the differential cross-sections for exponential centrality
# plt.figure(figsize=(8, 9))

plt.step(exp_centrality_bins_signal_0, dsigma_exp_centrality_signal_0, where="mid", alpha=0.7, label="LHeC@1.2 TeV : Signal ($w^+ w^-) [f_{M_0} / \Lambda^4$)", color="red", linewidth=3)
plt.step(exp_centrality_bins_background, dsigma_exp_centrality_background, where="mid", label="LHeC@1.2 TeV : SM background ($w^+ w^-$)", color="blue", linewidth=3)
plt.xlabel(r"$C_{\ell}^{\mathrm{exp}}$")
plt.ylabel(r"$\frac{d\sigma}{dC_{\ell}^{\mathrm{exp}}} \ \mathrm{[pb]}$")
plt.title(r"Delphes simulation : $e^- p \to e^- w^+ w^- p \to e^- j j \ell \nu_{\ell} p$ : LHeC@1.2 TeV", fontsize=20)
plt.yscale("log")
plt.legend()
plt.grid(True, linestyle="--", alpha=0.6)
# Add formula inside the plot
plt.text(0.5, 0.001, r"$C_{\ell}^{\mathrm{exp}} = e^{-|C_{\ell}|}$", color="black")
plt.tight_layout()
plt.savefig("/home/hamzeh-khanpour/Documents/GitHub/LHeC_Fast_Simulation/Pythia8_Delphes/differential_cross_section_exp_centrality.png", dpi=600)
plt.show()





# Plot the differential cross-sections for jet centrality
# plt.figure(figsize=(8, 9))

plt.step(jet_centrality_bins_signal_0, dsigma_jet_centrality_signal_0, where="mid", alpha=0.7, label="LHeC@1.2 TeV : Signal ($w^+ w^-) [f_{M_0} / \Lambda^4$)", color="red", linewidth=3)
plt.step(jet_centrality_bins_background, dsigma_jet_centrality_background, where="mid", alpha=0.7, label="LHeC@1.2 TeV : SM background ($w^+ w^-$)", color="blue", linewidth=3)
plt.xlabel(r"$C_{\mathrm{jets}}$")
plt.ylabel(r"$\frac{d\sigma}{dC_{\mathrm{jets}}} \ \mathrm{[pb]}$")
plt.title(r"Delphes simulation : $e^- p \to e^- w^+ w^- p \to e^- j j \ell \nu_{\ell} p$ : LHeC@1.2 TeV", fontsize=20)
plt.yscale("log")
plt.legend()
plt.grid(True, linestyle="--", alpha=0.6)
# Add formula inside the plot
plt.text(0.5, 0.001, r"$C_{\mathrm{jets}} = \frac{|\eta_{\mathrm{jet1}} + \eta_{\mathrm{jet2}}|}{2}$", color="black")
plt.tight_layout()
plt.savefig("/home/hamzeh-khanpour/Documents/GitHub/LHeC_Fast_Simulation/Pythia8_Delphes/differential_cross_section_jet_centrality.png", dpi=600)
plt.show()





# Plot the differential cross-sections for Δηjj
plt.step(delta_eta_jj_bins_signal_0, dsigma_delta_eta_jj_signal_0, where="mid", alpha=0.7, label="LHeC@1.2 TeV : Signal ($w^+ w^-) [f_{M_0} / \Lambda^4$)", color="red", linewidth=3)
plt.step(delta_eta_jj_bins_background, dsigma_delta_eta_jj_background, where="mid", alpha=0.7, label="LHeC@1.2 TeV : SM background ($w^+ w^-$)", color="blue", linewidth=3)
plt.xlabel(r"$\Delta \eta_{jj}$")
plt.ylabel(r"$\frac{d\sigma}{d\Delta \eta_{jj}} \ \mathrm{[pb]}$")
plt.title(r"Delphes simulation : $e^- p \to e^- w^+ w^- p \to e^- j j \ell \nu_{\ell} p$ : LHeC@1.2 TeV", fontsize=20)
plt.yscale("log")
plt.legend()
plt.grid(True, linestyle="--", alpha=0.6)
plt.tight_layout()
plt.savefig("/home/hamzeh-khanpour/Documents/GitHub/LHeC_Fast_Simulation/Pythia8_Delphes/differential_cross_section_delta_eta_jj.png", dpi=600)
plt.show()





#=========================================================================
#=========================================================================



# Open the output ROOT file
output_file = ROOT.TFile("/home/hamzeh-khanpour/Documents/GitHub/LHeC_Fast_Simulation/Pythia8_Delphes/output_histograms.root", "RECREATE")

# Create directories for signal and background
signal_dir = output_file.mkdir("signal_0")
background_dir = output_file.mkdir("SM_background")

# Save signal histograms to the "signal_0" branch
signal_dir.cd()  # Switch to the "signal_0" directory
hist_lepton_signal_0.Write()
hist_leading_jet_signal_0.Write()
hist_lepton_eta_signal_0.Write()
hist_delta_r_signal_0.Write()
hist_missing_et_signal_0.Write()
hist_centrality_signal_0.Write()
hist_exp_centrality_signal_0.Write()
hist_jet_centrality_signal_0.Write()
hist_delta_eta_jj_signal_0.Write()

# Save background histograms to the "SM_background" branch
background_dir.cd()  # Switch to the "SM_background" directory
hist_lepton_background.Write()
hist_leading_jet_background.Write()
hist_lepton_eta_background.Write()
hist_delta_r_background.Write()
hist_missing_et_background.Write()
hist_centrality_background.Write()
hist_exp_centrality_background.Write()
hist_jet_centrality_background.Write()
hist_delta_eta_jj_background.Write()

# Close the output file
output_file.Close()
print("Histograms saved to output_histograms.root with 'signal_0' and 'SM_background' branches.")



#=========================================================================
#=========================================================================



