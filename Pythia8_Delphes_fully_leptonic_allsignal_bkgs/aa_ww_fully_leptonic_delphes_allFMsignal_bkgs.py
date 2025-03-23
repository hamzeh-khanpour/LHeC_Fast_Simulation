

## ================================================================================
##        Hamzeh Khanpour  --- January 2025
## ================================================================================



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



# Path to the ROOT files : signal and background
signal_files = {
    "$FM_{0} / \Lambda^4$": "aa_ww_fully_leptonic_NP_1_FM0_100.root",
    "$FM_{1} / \Lambda^4$": "aa_ww_fully_leptonic_NP_1_FM1_100.root",
    "$FM_{2} / \Lambda^4$": "aa_ww_fully_leptonic_NP_1_FM2_100.root",
    "$FM_{3} / \Lambda^4$": "aa_ww_fully_leptonic_NP_1_FM3_100.root",
}

aa_ww_background_file_path     = "aa_ww_fully_leptonic_SM.root"
aa_ttbar_background_file_path  = "aa_ttbar_fully_leptonic.root"
aa_tautau_background_file_path = "aa_tautau_fully_leptonic.root"
aa_mumu_background_file_path   = "aa_mumu.root"



inclusive_ttbar_background_file_path   = "inclusive_ttbar_fully_leptonic.root"
single_top_background_file_path        = "single_top_fully_leptonic.root"
w_production_background_file_path      = "w_production_fully_leptonic.root"
z_production_background_file_path      = "z_production_fully_leptonic.root"



# Load Delphes library
ROOT.gSystem.Load("libDelphes")
ROOT.gInterpreter.Declare('#include "classes/DelphesClasses.h"')
ROOT.gInterpreter.Declare('#include "external/ExRootAnalysis/ExRootTreeReader.h"')


#=========================================================================
#=========================================================================


# Function to process a ROOT file and fill histograms for fully leptonic events
def process_file(
    file_path,
    hist_leading_lepton,
    hist_subleading_lepton,
    hist_m_ll,
    hist_pt_ll,
    hist_rapidity_ll,
    hist_lepton_separation,
    hist_azimuthal_angle_separation,
    hist_rapidity_difference,
    hist_met,
    hist_mT_W,  # ✅ Transverse Mass of W boson
    hist_m_WW,  # ✅ Reconstructed Diboson System Mass
    hist_pT_WW,  # ✅ Diboson Transverse Momentum
    hist_leading_lepton_eta,
    hist_subleading_lepton_eta,
):

    # Open the ROOT file
    chain = ROOT.TChain("Delphes")
    chain.Add(file_path)

    # Create ExRootTreeReader object
    treeReader = ROOT.ExRootTreeReader(chain)
    numberOfEntries =  10000   #   treeReader.GetEntries()

    # Counters for efficiency calculation
    total_events = numberOfEntries
    selected_events_pre = 0
    selected_events_final = 0

    # Get branches for electrons, muons, jets, and MET
    branchElectron = treeReader.UseBranch("Electron")
    branchMuon = treeReader.UseBranch("Muon")
    branchJet = treeReader.UseBranch("Jet")
    branchMissingET = treeReader.UseBranch("MissingET")


    # Process each event
    for entry in range(numberOfEntries):
        treeReader.ReadEntry(entry)


        # Collect leptons (electrons + muons)
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


        # Collect jets (for vetoing)
        jets = []
        for i in range(branchJet.GetEntries()):
            jet = branchJet.At(i)
            jet_vec = TLorentzVector()
            jet_vec.SetPtEtaPhiM(jet.PT, jet.Eta, jet.Phi, jet.Mass)
            if jet.PT > 10:
                jets.append(jet_vec)


        # Apply selection criteria: exactly two leptons and no jets
        if len(leptons) != 2 or len(jets) > 0:
            continue

        # Count selected events
        selected_events_pre += 1


        # Identify leading and subleading leptons
        leading_lepton = max(leptons, key=lambda lep: lep.Pt())
        subleading_lepton = min(leptons, key=lambda lep: lep.Pt())

        # Fill histograms for individual leptons
        hist_leading_lepton.Fill(leading_lepton.Pt())
        hist_subleading_lepton.Fill(subleading_lepton.Pt())
        hist_leading_lepton_eta.Fill(leading_lepton.Eta())
        hist_subleading_lepton_eta.Fill(subleading_lepton.Eta())

        # **Dilepton system calculations**
        dilepton_system = leading_lepton + subleading_lepton

        # Compute additional kinematic variables
        delta_r = leading_lepton.DeltaR(subleading_lepton)  # Lepton separation ΔR
        delta_phi = abs(leading_lepton.DeltaPhi(subleading_lepton))  # Azimuthal Angle Separation Δφ
        delta_y = abs(leading_lepton.Rapidity() - subleading_lepton.Rapidity())  # Rapidity Difference Δy


        # Fill histograms for dilepton kinematics
        hist_m_ll.Fill(dilepton_system.M())  # Invariant mass M(ll)
        hist_pt_ll.Fill(dilepton_system.Pt())  # Transverse momentum pT(ll)
        hist_rapidity_ll.Fill(dilepton_system.Rapidity())  # Rapidity Y(ll)


        # Fill histograms for dilepton kinematics
        hist_lepton_separation.Fill(delta_r)
        hist_azimuthal_angle_separation.Fill(delta_phi)
        hist_rapidity_difference.Fill(delta_y)

        # **Missing Transverse Energy (MET)**
        if branchMissingET.GetEntries() > 0:
            missing_et = branchMissingET.At(0)
            hist_met.Fill(missing_et.MET)  # Fill MET histogram

            # **Transverse Mass Calculation \( M_T(W) \)**
            M_T_W = np.sqrt(
                2 * leading_lepton.Pt() * missing_et.MET * (1 - np.cos(leading_lepton.Phi() - missing_et.Phi))
            )
            hist_mT_W.Fill(M_T_W)

            # **Diboson System Reconstruction**
            WW_system = dilepton_system + ROOT.TLorentzVector()
            WW_system.SetPxPyPzE(missing_et.MET * np.cos(missing_et.Phi),
                                  missing_et.MET * np.sin(missing_et.Phi),
                                  0,  # Approximating neutrino longitudinal momentum as unknown
                                  missing_et.MET)

            hist_m_WW.Fill(WW_system.M())  # Diboson system mass M(WW)
            hist_pT_WW.Fill(WW_system.Pt())  # Diboson transverse momentum pT(WW)


        # **✅ Final Event Selection **
        if leading_lepton.Pt() < 1.0  or  subleading_lepton.Pt() < 1.0:
            continue

        # Count selected events
        selected_events_final += 1


    # Calculate selection efficiency
    efficiency_pre = selected_events_pre / total_events if total_events > 0 else 0
    efficiency_final = selected_events_final / total_events if total_events > 0 else 0

    return (
        hist_leading_lepton,
        hist_subleading_lepton,
        hist_m_ll,
        hist_pt_ll,
        hist_rapidity_ll,
        hist_lepton_separation,
        hist_azimuthal_angle_separation,
        hist_rapidity_difference,
        hist_met,
        hist_mT_W,
        hist_m_WW,
        hist_pT_WW,
        hist_leading_lepton_eta,
        hist_subleading_lepton_eta,
        efficiency_pre,
        efficiency_final,
    )






#=========================================================================
#=========================================================================



# Parameters for differential cross-section

signal_cross_sections = {
    "$FM_{0} / \Lambda^4$": 0.00338968,   # pb
    "$FM_{1} / \Lambda^4$": 0.00343213,   # pb
    "$FM_{2} / \Lambda^4$": 0.00498621,   # pb
    "$FM_{3} / \Lambda^4$": 0.00376192    # pb
}

aa_ww_background_cross_section     = 0.00357101  # pb
aa_ttbar_background_cross_section  = 2.3846e-03  # pb  * 10^{+3}
aa_tautau_background_cross_section = 2.51510000  # pb
aa_mumu_background_cross_section   = 2.57270000  # pb


inclusive_ttbar_background_cross_section   = 0.00041134           # pb
single_top_background_cross_section        = 0.34049              # pb
w_production_background_cross_section      = 0.47776399999999997  # pb
z_production_background_cross_section      = 0.0764974266628      # pb




num_bins = 50

# ✅ Leading & Subleading Lepton Kinematics
pt_range_leading_lepton = (0, 200)     # Range for leading lepton pT
pt_range_subleading_lepton = (0, 200)  # Range for subleading lepton pT
eta_range_leading_lepton = (-5, 5)     # Range for leading lepton pseudorapidity
eta_range_subleading_lepton = (-5, 5)  # Range for subleading lepton pseudorapidity

# ✅ Dilepton System Kinematics
m_ll_range = (0, 200)         # Range for dilepton invariant mass M(ll)
pt_ll_range = (0, 200)         # Range for dilepton transverse momentum pT(ll)
rapidity_ll_range = (-5, 5)    # Range for dilepton system rapidity Y(ll)

# ✅ Angular Separations
delta_r_range = (0, 6)         # Range for ΔR between leptons
delta_phi_range = (0, np.pi)   # Range for Δφ between leptons
delta_y_range = (0, 5)         # Range for rapidity difference Δy between leptons

# ✅ Missing Energy & Neutrino Variables
met_range = (0, 200)           # Range for Missing Transverse Energy (MET)

# ✅ W Boson & Diboson System
mT_W_range = (0, 200)          # Range for transverse mass of W boson
m_WW_range = (0, 200)         # Range for diboson system mass M(WW)
pT_WW_range = (0, 200)         # Range for transverse momentum of diboson system pT(WW)



# ✅ Calculate Bin Widths
bin_width_pt_leading_lepton = (pt_range_leading_lepton[1] - pt_range_leading_lepton[0]) / num_bins
bin_width_pt_subleading_lepton = (pt_range_subleading_lepton[1] - pt_range_subleading_lepton[0]) / num_bins
bin_width_eta_leading_lepton = (eta_range_leading_lepton[1] - eta_range_leading_lepton[0]) / num_bins
bin_width_eta_subleading_lepton = (eta_range_subleading_lepton[1] - eta_range_subleading_lepton[0]) / num_bins

bin_width_m_ll = (m_ll_range[1] - m_ll_range[0]) / num_bins
bin_width_pt_ll = (pt_ll_range[1] - pt_ll_range[0]) / num_bins
bin_width_rapidity_ll = (rapidity_ll_range[1] - rapidity_ll_range[0]) / num_bins

bin_width_delta_r = (delta_r_range[1] - delta_r_range[0]) / num_bins
bin_width_delta_phi = (delta_phi_range[1] - delta_phi_range[0]) / num_bins
bin_width_delta_y = (delta_y_range[1] - delta_y_range[0]) / num_bins

bin_width_met = (met_range[1] - met_range[0]) / num_bins

bin_width_mT_W = (mT_W_range[1] - mT_W_range[0]) / num_bins
bin_width_m_WW = (m_WW_range[1] - m_WW_range[0]) / num_bins
bin_width_pT_WW = (pT_WW_range[1] - pT_WW_range[0]) / num_bins






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




"""

# Function to calculate differential cross-section from a ROOT histogram
def calculate_dsigma(histogram, total_cross_section, bin_width):

    #Converts a ROOT TH1 histogram into a normalized differential cross-section.

    #Args:
        #histogram (ROOT.TH1F): The histogram to process.
        #total_cross_section (float): The total cross-section for this process (in pb).
        #bin_width (float): The fixed bin width for normalization (in GeV or appropriate units).

    #Returns:
        #Tuple of (bin_centers, dsigma_values)


    # Safety check for empty histogram
    if histogram.GetEntries() == 0:
        print(f"⚠️  Warning: Histogram '{histogram.GetName()}' is empty. Returning zeros.")
        num_bins = histogram.GetNbinsX()
        bin_centers = [histogram.GetBinCenter(i) for i in range(1, num_bins + 1)]
        dsigma = [0.0] * num_bins
        return bin_centers, dsigma

    # Extract counts and compute differential cross-section
    num_bins = histogram.GetNbinsX()
    total_counts = histogram.Integral()
    bin_centers = [histogram.GetBinCenter(i) for i in range(1, num_bins + 1)]
    dsigma = [
        histogram.GetBinContent(i) * (total_cross_section / total_counts) / bin_width
        for i in range(1, num_bins + 1)
    ]

    return bin_centers, dsigma

"""


#=========================================================================
#=========================================================================

# Dictionary to store histograms for each signal
signal_histograms = {}

for signal_name in signal_files:
    signal_histograms[signal_name] = {
        # ✅ Leading & Subleading Lepton Kinematics
        "hist_leading_lepton": ROOT.TH1F(f"hist_leading_lepton_{signal_name}", f"Leading Lepton pT ({signal_name}); p_{{T}}^{{leading}} [GeV]; Entries", num_bins, *pt_range_leading_lepton),
        "hist_subleading_lepton": ROOT.TH1F(f"hist_subleading_lepton_{signal_name}", f"Subleading Lepton pT ({signal_name}); p_{{T}}^{{subleading}} [GeV]; Entries", num_bins, *pt_range_subleading_lepton),
        "hist_leading_lepton_eta": ROOT.TH1F(f"hist_leading_lepton_eta_{signal_name}", f"Leading Lepton Eta ({signal_name}); #eta_{{leading}}; Entries", num_bins, *eta_range_leading_lepton),
        "hist_subleading_lepton_eta": ROOT.TH1F(f"hist_subleading_lepton_eta_{signal_name}", f"Subleading Lepton Eta ({signal_name}); #eta_{{subleading}}; Entries", num_bins, *eta_range_subleading_lepton),

        # ✅ Dilepton System Kinematics
        "hist_m_ll": ROOT.TH1F(f"hist_m_ll_{signal_name}", f"Dilepton Invariant Mass ({signal_name}); M_{{ll}} [GeV]; Entries", num_bins, *m_ll_range),
        "hist_pt_ll": ROOT.TH1F(f"hist_pt_ll_{signal_name}", f"Dilepton pT ({signal_name}); p_{{T}}^{{ll}} [GeV]; Entries", num_bins, *pt_ll_range),
        "hist_rapidity_ll": ROOT.TH1F(f"hist_rapidity_ll_{signal_name}", f"Dilepton Rapidity ({signal_name}); Y_{{ll}}; Entries", num_bins, *rapidity_ll_range),

        # ✅ Angular Separations
        "hist_lepton_separation": ROOT.TH1F(f"hist_lepton_separation_{signal_name}", f"Lepton Separation ΔR ({signal_name}); ΔR_{{ll}}; Entries", num_bins, *delta_r_range),
        "hist_azimuthal_angle_separation": ROOT.TH1F(f"hist_azimuthal_angle_separation_{signal_name}", f"Azimuthal Angle Separation ({signal_name}); Δφ_{{ll}}; Entries", num_bins, *delta_phi_range),
        "hist_rapidity_difference": ROOT.TH1F(f"hist_rapidity_difference_{signal_name}", f"Rapidity Difference Δy ({signal_name}); Δy_{{ll}}; Entries", num_bins, *delta_y_range),

        # ✅ Missing Transverse Energy
        "hist_met": ROOT.TH1F(f"hist_met_{signal_name}", f"Missing ET ({signal_name}); MET [GeV]; Entries", num_bins, *met_range),

        # ✅ W Boson & Diboson System
        "hist_mT_W": ROOT.TH1F(f"hist_mT_W_{signal_name}", f"Transverse Mass of W ({signal_name}); M_T^{{W}} [GeV]; Entries", num_bins, *mT_W_range),
        "hist_m_WW": ROOT.TH1F(f"hist_m_WW_{signal_name}", f"Diboson Mass ({signal_name}); M_{{WW}} [GeV]; Entries", num_bins, *m_WW_range),
        "hist_pT_WW": ROOT.TH1F(f"hist_pT_WW_{signal_name}", f"Diboson pT ({signal_name}); p_{{T}}^{{WW}} [GeV]; Entries", num_bins, *pT_WW_range),
    }



# ✅ SM WW Background histograms for fully leptonic analysis
hist_leading_lepton_background = ROOT.TH1F("hist_leading_lepton_background", "Leading Lepton pT (Background); p_{T}^{leading} [GeV]; Entries", num_bins, *pt_range_leading_lepton)
hist_subleading_lepton_background = ROOT.TH1F("hist_subleading_lepton_background", "Subleading Lepton pT (Background); p_{T}^{subleading} [GeV]; Entries", num_bins, *pt_range_subleading_lepton)
hist_leading_lepton_eta_background = ROOT.TH1F("hist_leading_lepton_eta_background", "Leading Lepton Eta (Background); #eta_{leading}; Entries", num_bins, *eta_range_leading_lepton)
hist_subleading_lepton_eta_background = ROOT.TH1F("hist_subleading_lepton_eta_background", "Subleading Lepton Eta (Background); #eta_{subleading}; Entries", num_bins, *eta_range_subleading_lepton)

hist_m_ll_background = ROOT.TH1F("hist_m_ll_background", "Dilepton Invariant Mass (Background); M_{ll} [GeV]; Entries", num_bins, *m_ll_range)
hist_pt_ll_background = ROOT.TH1F("hist_pt_ll_background", "Dilepton pT (Background); p_{T}^{ll} [GeV]; Entries", num_bins, *pt_ll_range)
hist_rapidity_ll_background = ROOT.TH1F("hist_rapidity_ll_background", "Dilepton Rapidity (Background); Y_{ll}; Entries", num_bins, *rapidity_ll_range)

hist_lepton_separation_background = ROOT.TH1F("hist_lepton_separation_background", "Lepton Separation ΔR (Background); ΔR_{ll}; Entries", num_bins, *delta_r_range)
hist_azimuthal_angle_separation_background = ROOT.TH1F("hist_azimuthal_angle_separation_background", "Azimuthal Angle Separation (Background); Δφ_{ll}; Entries", num_bins, *delta_phi_range)
hist_rapidity_difference_background = ROOT.TH1F("hist_rapidity_difference_background", "Rapidity Difference Δy (Background); Δy_{ll}; Entries", num_bins, *delta_y_range)

hist_met_background = ROOT.TH1F("hist_met_background", "Missing ET (Background); MET [GeV]; Entries", num_bins, *met_range)

hist_mT_W_background = ROOT.TH1F("hist_mT_W_background", "Transverse Mass of W (Background); M_T^{W} [GeV]; Entries", num_bins, *mT_W_range)
hist_m_WW_background = ROOT.TH1F("hist_m_WW_background", "Diboson Mass (Background); M_{WW} [GeV]; Entries", num_bins, *m_WW_range)
hist_pT_WW_background = ROOT.TH1F("hist_pT_WW_background", "Diboson pT (Background); p_{T}^{WW} [GeV]; Entries", num_bins, *pT_WW_range)




# ✅ ttbar Background Histograms
hist_leading_lepton_ttbar = ROOT.TH1F("hist_leading_lepton_ttbar", "Leading Lepton pT (ttbar); p_{T}^{leading} [GeV]; Entries", num_bins, *pt_range_leading_lepton)
hist_subleading_lepton_ttbar = ROOT.TH1F("hist_subleading_lepton_ttbar", "Subleading Lepton pT (ttbar); p_{T}^{subleading} [GeV]; Entries", num_bins, *pt_range_subleading_lepton)
hist_leading_lepton_eta_ttbar = ROOT.TH1F("hist_leading_lepton_eta_ttbar", "Leading Lepton Eta (ttbar); #eta_{leading}; Entries", num_bins, *eta_range_leading_lepton)
hist_subleading_lepton_eta_ttbar = ROOT.TH1F("hist_subleading_lepton_eta_ttbar", "Subleading Lepton Eta (ttbar); #eta_{subleading}; Entries", num_bins, *eta_range_subleading_lepton)

hist_m_ll_ttbar = ROOT.TH1F("hist_m_ll_ttbar", "Dilepton Invariant Mass (ttbar); M_{ll} [GeV]; Entries", num_bins, *m_ll_range)
hist_pt_ll_ttbar = ROOT.TH1F("hist_pt_ll_ttbar", "Dilepton pT (ttbar); p_{T}^{ll} [GeV]; Entries", num_bins, *pt_ll_range)
hist_rapidity_ll_ttbar = ROOT.TH1F("hist_rapidity_ll_ttbar", "Dilepton Rapidity (ttbar); Y_{ll}; Entries", num_bins, *rapidity_ll_range)

hist_lepton_separation_ttbar = ROOT.TH1F("hist_lepton_separation_ttbar", "Lepton Separation ΔR (ttbar); ΔR_{ll}; Entries", num_bins, *delta_r_range)
hist_azimuthal_angle_separation_ttbar = ROOT.TH1F("hist_azimuthal_angle_separation_ttbar", "Azimuthal Angle Separation (ttbar); Δφ_{ll}; Entries", num_bins, *delta_phi_range)
hist_rapidity_difference_ttbar = ROOT.TH1F("hist_rapidity_difference_ttbar", "Rapidity Difference Δy (ttbar); Δy_{ll}; Entries", num_bins, *delta_y_range)

hist_met_ttbar = ROOT.TH1F("hist_met_ttbar", "Missing ET (ttbar); MET [GeV]; Entries", num_bins, *met_range)
hist_mT_W_ttbar = ROOT.TH1F("hist_mT_W_ttbar", "Transverse Mass of W (ttbar); M_T^{W} [GeV]; Entries", num_bins, *mT_W_range)
hist_m_WW_ttbar = ROOT.TH1F("hist_m_WW_ttbar", "Diboson Mass (ttbar); M_{WW} [GeV]; Entries", num_bins, *m_WW_range)
hist_pT_WW_ttbar = ROOT.TH1F("hist_pT_WW_ttbar", "Diboson pT (ttbar); p_{T}^{WW} [GeV]; Entries", num_bins, *pT_WW_range)






# ✅ tautau Background Histograms
hist_leading_lepton_tautau = ROOT.TH1F("hist_leading_lepton_tautau", "Leading Lepton pT (τ⁺τ⁻); p_{T}^{leading} [GeV]; Entries", num_bins, *pt_range_leading_lepton)
hist_subleading_lepton_tautau = ROOT.TH1F("hist_subleading_lepton_tautau", "Subleading Lepton pT (τ⁺τ⁻); p_{T}^{subleading} [GeV]; Entries", num_bins, *pt_range_subleading_lepton)
hist_leading_lepton_eta_tautau = ROOT.TH1F("hist_leading_lepton_eta_tautau", "Leading Lepton Eta (τ⁺τ⁻); #eta_{leading}; Entries", num_bins, *eta_range_leading_lepton)
hist_subleading_lepton_eta_tautau = ROOT.TH1F("hist_subleading_lepton_eta_tautau", "Subleading Lepton Eta (τ⁺τ⁻); #eta_{subleading}; Entries", num_bins, *eta_range_subleading_lepton)

hist_m_ll_tautau = ROOT.TH1F("hist_m_ll_tautau", "Dilepton Invariant Mass (τ⁺τ⁻); M_{ll} [GeV]; Entries", num_bins, *m_ll_range)
hist_pt_ll_tautau = ROOT.TH1F("hist_pt_ll_tautau", "Dilepton pT (τ⁺τ⁻); p_{T}^{ll} [GeV]; Entries", num_bins, *pt_ll_range)
hist_rapidity_ll_tautau = ROOT.TH1F("hist_rapidity_ll_tautau", "Dilepton Rapidity (τ⁺τ⁻); Y_{ll}; Entries", num_bins, *rapidity_ll_range)

hist_lepton_separation_tautau = ROOT.TH1F("hist_lepton_separation_tautau", "Lepton Separation ΔR (τ⁺τ⁻); ΔR_{ll}; Entries", num_bins, *delta_r_range)
hist_azimuthal_angle_separation_tautau = ROOT.TH1F("hist_azimuthal_angle_separation_tautau", "Azimuthal Angle Separation (τ⁺τ⁻); Δφ_{ll}; Entries", num_bins, *delta_phi_range)
hist_rapidity_difference_tautau = ROOT.TH1F("hist_rapidity_difference_tautau", "Rapidity Difference Δy (τ⁺τ⁻); Δy_{ll}; Entries", num_bins, *delta_y_range)

hist_met_tautau = ROOT.TH1F("hist_met_tautau", "Missing ET (τ⁺τ⁻); MET [GeV]; Entries", num_bins, *met_range)
hist_mT_W_tautau = ROOT.TH1F("hist_mT_W_tautau", "Transverse Mass of W (τ⁺τ⁻); M_T^{W} [GeV]; Entries", num_bins, *mT_W_range)
hist_m_WW_tautau = ROOT.TH1F("hist_m_WW_tautau", "Diboson Mass (τ⁺τ⁻); M_{WW} [GeV]; Entries", num_bins, *m_WW_range)
hist_pT_WW_tautau = ROOT.TH1F("hist_pT_WW_tautau", "Diboson pT (τ⁺τ⁻); p_{T}^{WW} [GeV]; Entries", num_bins, *pT_WW_range)






# ✅ mumu Background Histograms
hist_leading_lepton_mumu = ROOT.TH1F("hist_leading_lepton_mumu", "Leading Lepton pT (μ⁺μ⁻); p_{T}^{leading} [GeV]; Entries", num_bins, *pt_range_leading_lepton)
hist_subleading_lepton_mumu = ROOT.TH1F("hist_subleading_lepton_mumu", "Subleading Lepton pT (μ⁺μ⁻); p_{T}^{subleading} [GeV]; Entries", num_bins, *pt_range_subleading_lepton)
hist_leading_lepton_eta_mumu = ROOT.TH1F("hist_leading_lepton_eta_mumu", "Leading Lepton Eta (μ⁺μ⁻); #eta_{leading}; Entries", num_bins, *eta_range_leading_lepton)
hist_subleading_lepton_eta_mumu = ROOT.TH1F("hist_subleading_lepton_eta_mumu", "Subleading Lepton Eta (μ⁺μ⁻); #eta_{subleading}; Entries", num_bins, *eta_range_subleading_lepton)

hist_m_ll_mumu = ROOT.TH1F("hist_m_ll_mumu", "Dilepton Invariant Mass (μ⁺μ⁻); M_{ll} [GeV]; Entries", num_bins, *m_ll_range)
hist_pt_ll_mumu = ROOT.TH1F("hist_pt_ll_mumu", "Dilepton pT (μ⁺μ⁻); p_{T}^{ll} [GeV]; Entries", num_bins, *pt_ll_range)
hist_rapidity_ll_mumu = ROOT.TH1F("hist_rapidity_ll_mumu", "Dilepton Rapidity (μ⁺μ⁻); Y_{ll}; Entries", num_bins, *rapidity_ll_range)

hist_lepton_separation_mumu = ROOT.TH1F("hist_lepton_separation_mumu", "Lepton Separation ΔR (μ⁺μ⁻); ΔR_{ll}; Entries", num_bins, *delta_r_range)
hist_azimuthal_angle_separation_mumu = ROOT.TH1F("hist_azimuthal_angle_separation_mumu", "Azimuthal Angle Separation (μ⁺μ⁻); Δφ_{ll}; Entries", num_bins, *delta_phi_range)
hist_rapidity_difference_mumu = ROOT.TH1F("hist_rapidity_difference_mumu", "Rapidity Difference Δy (μ⁺μ⁻); Δy_{ll}; Entries", num_bins, *delta_y_range)

hist_met_mumu = ROOT.TH1F("hist_met_mumu", "Missing ET (μ⁺μ⁻); MET [GeV]; Entries", num_bins, *met_range)
hist_mT_W_mumu = ROOT.TH1F("hist_mT_W_mumu", "Transverse Mass of W (μ⁺μ⁻); M_T^{W} [GeV]; Entries", num_bins, *mT_W_range)
hist_m_WW_mumu = ROOT.TH1F("hist_m_WW_mumu", "Diboson Mass (μ⁺μ⁻); M_{WW} [GeV]; Entries", num_bins, *m_WW_range)
hist_pT_WW_mumu = ROOT.TH1F("hist_pT_WW_mumu", "Diboson pT (μ⁺μ⁻); p_{T}^{WW} [GeV]; Entries", num_bins, *pT_WW_range)




# ✅ inclusive ttbar Background Histograms

hist_leading_lepton_ttbar_inc = ROOT.TH1F("hist_leading_lepton_ttbar_inc", "Leading Lepton pT (inclusive ttbar)", num_bins, *pt_range_leading_lepton)
hist_subleading_lepton_ttbar_inc = ROOT.TH1F("hist_subleading_lepton_ttbar_inc", "Subleading Lepton pT (inclusive ttbar)", num_bins, *pt_range_subleading_lepton)
hist_leading_lepton_eta_ttbar_inc = ROOT.TH1F("hist_leading_lepton_eta_ttbar_inc", "Leading Lepton Eta (inclusive ttbar)", num_bins, *eta_range_leading_lepton)
hist_subleading_lepton_eta_ttbar_inc = ROOT.TH1F("hist_subleading_lepton_eta_ttbar_inc", "Subleading Lepton Eta (inclusive ttbar)", num_bins, *eta_range_subleading_lepton)

hist_m_ll_ttbar_inc = ROOT.TH1F("hist_m_ll_ttbar_inc", "Dilepton Invariant Mass (inclusive ttbar)", num_bins, *m_ll_range)
hist_pt_ll_ttbar_inc = ROOT.TH1F("hist_pt_ll_ttbar_inc", "Dilepton pT (inclusive ttbar)", num_bins, *pt_ll_range)
hist_rapidity_ll_ttbar_inc = ROOT.TH1F("hist_rapidity_ll_ttbar_inc", "Dilepton Rapidity (inclusive ttbar)", num_bins, *rapidity_ll_range)

hist_lepton_separation_ttbar_inc = ROOT.TH1F("hist_lepton_separation_ttbar_inc", "Lepton Separation ΔR (inclusive ttbar)", num_bins, *delta_r_range)
hist_azimuthal_angle_separation_ttbar_inc = ROOT.TH1F("hist_azimuthal_angle_separation_ttbar_inc", "Azimuthal Angle Separation (inclusive ttbar)", num_bins, *delta_phi_range)
hist_rapidity_difference_ttbar_inc = ROOT.TH1F("hist_rapidity_difference_ttbar_inc", "Rapidity Difference Δy (inclusive ttbar)", num_bins, *delta_y_range)

hist_met_ttbar_inc = ROOT.TH1F("hist_met_ttbar_inc", "Missing ET (inclusive ttbar)", num_bins, *met_range)
hist_mT_W_ttbar_inc = ROOT.TH1F("hist_mT_W_ttbar_inc", "Transverse Mass of W (inclusive ttbar)", num_bins, *mT_W_range)
hist_m_WW_ttbar_inc = ROOT.TH1F("hist_m_WW_ttbar_inc", "Diboson Mass (inclusive ttbar)", num_bins, *m_WW_range)
hist_pT_WW_ttbar_inc = ROOT.TH1F("hist_pT_WW_ttbar_inc", "Diboson pT (inclusive ttbar)", num_bins, *pT_WW_range)





# ✅  Single Top Background Histograms

hist_leading_lepton_single_top = ROOT.TH1F("hist_leading_lepton_single_top", "Leading Lepton pT (single top)", num_bins, *pt_range_leading_lepton)
hist_subleading_lepton_single_top = ROOT.TH1F("hist_subleading_lepton_single_top", "Subleading Lepton pT (single top)", num_bins, *pt_range_subleading_lepton)
hist_leading_lepton_eta_single_top = ROOT.TH1F("hist_leading_lepton_eta_single_top", "Leading Lepton Eta (single top)", num_bins, *eta_range_leading_lepton)
hist_subleading_lepton_eta_single_top = ROOT.TH1F("hist_subleading_lepton_eta_single_top", "Subleading Lepton Eta (single top)", num_bins, *eta_range_subleading_lepton)

hist_m_ll_single_top = ROOT.TH1F("hist_m_ll_single_top", "Dilepton Invariant Mass (single top)", num_bins, *m_ll_range)
hist_pt_ll_single_top = ROOT.TH1F("hist_pt_ll_single_top", "Dilepton pT (single top)", num_bins, *pt_ll_range)
hist_rapidity_ll_single_top = ROOT.TH1F("hist_rapidity_ll_single_top", "Dilepton Rapidity (single top)", num_bins, *rapidity_ll_range)

hist_lepton_separation_single_top = ROOT.TH1F("hist_lepton_separation_single_top", "Lepton Separation ΔR (single top)", num_bins, *delta_r_range)
hist_azimuthal_angle_separation_single_top = ROOT.TH1F("hist_azimuthal_angle_separation_single_top", "Azimuthal Angle Separation (single top)", num_bins, *delta_phi_range)
hist_rapidity_difference_single_top = ROOT.TH1F("hist_rapidity_difference_single_top", "Rapidity Difference Δy (single top)", num_bins, *delta_y_range)

hist_met_single_top = ROOT.TH1F("hist_met_single_top", "Missing ET (single top)", num_bins, *met_range)
hist_mT_W_single_top = ROOT.TH1F("hist_mT_W_single_top", "Transverse Mass of W (single top)", num_bins, *mT_W_range)
hist_m_WW_single_top = ROOT.TH1F("hist_m_WW_single_top", "Diboson Mass (single top)", num_bins, *m_WW_range)
hist_pT_WW_single_top = ROOT.TH1F("hist_pT_WW_single_top", "Diboson pT (single top)", num_bins, *pT_WW_range)




# ✅   W Production Background Histograms

hist_leading_lepton_w = ROOT.TH1F("hist_leading_lepton_w", "Leading Lepton pT (W production)", num_bins, *pt_range_leading_lepton)
hist_subleading_lepton_w = ROOT.TH1F("hist_subleading_lepton_w", "Subleading Lepton pT (W production)", num_bins, *pt_range_subleading_lepton)
hist_leading_lepton_eta_w = ROOT.TH1F("hist_leading_lepton_eta_w", "Leading Lepton Eta (W production)", num_bins, *eta_range_leading_lepton)
hist_subleading_lepton_eta_w = ROOT.TH1F("hist_subleading_lepton_eta_w", "Subleading Lepton Eta (W production)", num_bins, *eta_range_subleading_lepton)

hist_m_ll_w = ROOT.TH1F("hist_m_ll_w", "Dilepton Invariant Mass (W production)", num_bins, *m_ll_range)
hist_pt_ll_w = ROOT.TH1F("hist_pt_ll_w", "Dilepton pT (W production)", num_bins, *pt_ll_range)
hist_rapidity_ll_w = ROOT.TH1F("hist_rapidity_ll_w", "Dilepton Rapidity (W production)", num_bins, *rapidity_ll_range)

hist_lepton_separation_w = ROOT.TH1F("hist_lepton_separation_w", "Lepton Separation ΔR (W production)", num_bins, *delta_r_range)
hist_azimuthal_angle_separation_w = ROOT.TH1F("hist_azimuthal_angle_separation_w", "Azimuthal Angle Separation (W production)", num_bins, *delta_phi_range)
hist_rapidity_difference_w = ROOT.TH1F("hist_rapidity_difference_w", "Rapidity Difference Δy (W production)", num_bins, *delta_y_range)

hist_met_w = ROOT.TH1F("hist_met_w", "Missing ET (W production)", num_bins, *met_range)
hist_mT_W_w = ROOT.TH1F("hist_mT_W_w", "Transverse Mass of W (W production)", num_bins, *mT_W_range)
hist_m_WW_w = ROOT.TH1F("hist_m_WW_w", "Diboson Mass (W production)", num_bins, *m_WW_range)
hist_pT_WW_w = ROOT.TH1F("hist_pT_WW_w", "Diboson pT (W production)", num_bins, *pT_WW_range)




# ✅   Z Production Background Histograms

hist_leading_lepton_z = ROOT.TH1F("hist_leading_lepton_z", "Leading Lepton pT (Z production)", num_bins, *pt_range_leading_lepton)
hist_subleading_lepton_z = ROOT.TH1F("hist_subleading_lepton_z", "Subleading Lepton pT (Z production)", num_bins, *pt_range_subleading_lepton)
hist_leading_lepton_eta_z = ROOT.TH1F("hist_leading_lepton_eta_z", "Leading Lepton Eta (Z production)", num_bins, *eta_range_leading_lepton)
hist_subleading_lepton_eta_z = ROOT.TH1F("hist_subleading_lepton_eta_z", "Subleading Lepton Eta (Z production)", num_bins, *eta_range_subleading_lepton)

hist_m_ll_z = ROOT.TH1F("hist_m_ll_z", "Dilepton Invariant Mass (Z production)", num_bins, *m_ll_range)
hist_pt_ll_z = ROOT.TH1F("hist_pt_ll_z", "Dilepton pT (Z production)", num_bins, *pt_ll_range)
hist_rapidity_ll_z = ROOT.TH1F("hist_rapidity_ll_z", "Dilepton Rapidity (Z production)", num_bins, *rapidity_ll_range)

hist_lepton_separation_z = ROOT.TH1F("hist_lepton_separation_z", "Lepton Separation ΔR (Z production)", num_bins, *delta_r_range)
hist_azimuthal_angle_separation_z = ROOT.TH1F("hist_azimuthal_angle_separation_z", "Azimuthal Angle Separation (Z production)", num_bins, *delta_phi_range)
hist_rapidity_difference_z = ROOT.TH1F("hist_rapidity_difference_z", "Rapidity Difference Δy (Z production)", num_bins, *delta_y_range)

hist_met_z = ROOT.TH1F("hist_met_z", "Missing ET (Z production)", num_bins, *met_range)
hist_mT_W_z = ROOT.TH1F("hist_mT_W_z", "Transverse Mass of W (Z production)", num_bins, *mT_W_range)
hist_m_WW_z = ROOT.TH1F("hist_m_WW_z", "Diboson Mass (Z production)", num_bins, *m_WW_range)
hist_pT_WW_z = ROOT.TH1F("hist_pT_WW_z", "Diboson pT (Z production)", num_bins, *pT_WW_range)








# Dictionary to store efficiencies for each signal
signal_efficiencies = {}

# Process all signal files dynamically
for signal_name, file_path in signal_files.items():
    print(f"Processing signal: {signal_name}")

    # Get the corresponding histograms for this signal
    histograms = signal_histograms[signal_name]



    # Process the signal file
    (histograms["hist_leading_lepton"], histograms["hist_subleading_lepton"], histograms["hist_m_ll"],
     histograms["hist_pt_ll"], histograms["hist_rapidity_ll"], histograms["hist_lepton_separation"],
     histograms["hist_azimuthal_angle_separation"], histograms["hist_rapidity_difference"], histograms["hist_met"],
     histograms["hist_mT_W"], histograms["hist_m_WW"], histograms["hist_pT_WW"], histograms["hist_leading_lepton_eta"],
     histograms["hist_subleading_lepton_eta"], efficiency_pre, efficiency_final) = process_file(
        file_path, histograms["hist_leading_lepton"], histograms["hist_subleading_lepton"], histograms["hist_m_ll"],
        histograms["hist_pt_ll"], histograms["hist_rapidity_ll"], histograms["hist_lepton_separation"],
        histograms["hist_azimuthal_angle_separation"], histograms["hist_rapidity_difference"], histograms["hist_met"],
        histograms["hist_mT_W"], histograms["hist_m_WW"], histograms["hist_pT_WW"], histograms["hist_leading_lepton_eta"],
        histograms["hist_subleading_lepton_eta"]
    )



    # Store efficiencies for this signal
    signal_efficiencies[signal_name] = {
        "efficiency_pre": efficiency_pre,
        "efficiency_final": efficiency_final
    }



# Process the SM WW background file separately
(hist_leading_lepton_background, hist_subleading_lepton_background, hist_m_ll_background,
 hist_pt_ll_background, hist_rapidity_ll_background, hist_lepton_separation_background,
 hist_azimuthal_angle_separation_background, hist_rapidity_difference_background, hist_met_background,
 hist_mT_W_background, hist_m_WW_background, hist_pT_WW_background, hist_leading_lepton_eta_background,
 hist_subleading_lepton_eta_background, background_efficiency_pre, background_efficiency_final) = process_file(
    aa_ww_background_file_path, hist_leading_lepton_background, hist_subleading_lepton_background, hist_m_ll_background,
    hist_pt_ll_background, hist_rapidity_ll_background, hist_lepton_separation_background,
    hist_azimuthal_angle_separation_background, hist_rapidity_difference_background, hist_met_background,
    hist_mT_W_background, hist_m_WW_background, hist_pT_WW_background, hist_leading_lepton_eta_background,
    hist_subleading_lepton_eta_background
)




(hist_leading_lepton_ttbar, hist_subleading_lepton_ttbar, hist_m_ll_ttbar,
 hist_pt_ll_ttbar, hist_rapidity_ll_ttbar, hist_lepton_separation_ttbar,
 hist_azimuthal_angle_separation_ttbar, hist_rapidity_difference_ttbar, hist_met_ttbar,
 hist_mT_W_ttbar, hist_m_WW_ttbar, hist_pT_WW_ttbar, hist_leading_lepton_eta_ttbar,
 hist_subleading_lepton_eta_ttbar, ttbar_efficiency_pre, ttbar_efficiency_final) = process_file(
    aa_ttbar_background_file_path, hist_leading_lepton_ttbar, hist_subleading_lepton_ttbar, hist_m_ll_ttbar,
    hist_pt_ll_ttbar, hist_rapidity_ll_ttbar, hist_lepton_separation_ttbar,
    hist_azimuthal_angle_separation_ttbar, hist_rapidity_difference_ttbar, hist_met_ttbar,
    hist_mT_W_ttbar, hist_m_WW_ttbar, hist_pT_WW_ttbar, hist_leading_lepton_eta_ttbar,
    hist_subleading_lepton_eta_ttbar
)




(hist_leading_lepton_tautau, hist_subleading_lepton_tautau, hist_m_ll_tautau,
 hist_pt_ll_tautau, hist_rapidity_ll_tautau, hist_lepton_separation_tautau,
 hist_azimuthal_angle_separation_tautau, hist_rapidity_difference_tautau, hist_met_tautau,
 hist_mT_W_tautau, hist_m_WW_tautau, hist_pT_WW_tautau, hist_leading_lepton_eta_tautau,
 hist_subleading_lepton_eta_tautau, tautau_efficiency_pre, tautau_efficiency_final) = process_file(
    aa_tautau_background_file_path, hist_leading_lepton_tautau, hist_subleading_lepton_tautau, hist_m_ll_tautau,
    hist_pt_ll_tautau, hist_rapidity_ll_tautau, hist_lepton_separation_tautau,
    hist_azimuthal_angle_separation_tautau, hist_rapidity_difference_tautau, hist_met_tautau,
    hist_mT_W_tautau, hist_m_WW_tautau, hist_pT_WW_tautau, hist_leading_lepton_eta_tautau,
    hist_subleading_lepton_eta_tautau
)





(hist_leading_lepton_mumu, hist_subleading_lepton_mumu, hist_m_ll_mumu,
 hist_pt_ll_mumu, hist_rapidity_ll_mumu, hist_lepton_separation_mumu,
 hist_azimuthal_angle_separation_mumu, hist_rapidity_difference_mumu, hist_met_mumu,
 hist_mT_W_mumu, hist_m_WW_mumu, hist_pT_WW_mumu, hist_leading_lepton_eta_mumu,
 hist_subleading_lepton_eta_mumu, mumu_efficiency_pre, mumu_efficiency_final) = process_file(
    aa_mumu_background_file_path, hist_leading_lepton_mumu, hist_subleading_lepton_mumu, hist_m_ll_mumu,
    hist_pt_ll_mumu, hist_rapidity_ll_mumu, hist_lepton_separation_mumu,
    hist_azimuthal_angle_separation_mumu, hist_rapidity_difference_mumu, hist_met_mumu,
    hist_mT_W_mumu, hist_m_WW_mumu, hist_pT_WW_mumu, hist_leading_lepton_eta_mumu,
    hist_subleading_lepton_eta_mumu
)




(hist_leading_lepton_ttbar_inc, hist_subleading_lepton_ttbar_inc, hist_m_ll_ttbar_inc,
 hist_pt_ll_ttbar_inc, hist_rapidity_ll_ttbar_inc, hist_lepton_separation_ttbar_inc,
 hist_azimuthal_angle_separation_ttbar_inc, hist_rapidity_difference_ttbar_inc, hist_met_ttbar_inc,
 hist_mT_W_ttbar_inc, hist_m_WW_ttbar_inc, hist_pT_WW_ttbar_inc, hist_leading_lepton_eta_ttbar_inc,
 hist_subleading_lepton_eta_ttbar_inc, ttbar_inc_efficiency_pre, ttbar_inc_efficiency_final) = process_file(
    inclusive_ttbar_background_file_path, hist_leading_lepton_ttbar_inc, hist_subleading_lepton_ttbar_inc, hist_m_ll_ttbar_inc,
    hist_pt_ll_ttbar_inc, hist_rapidity_ll_ttbar_inc, hist_lepton_separation_ttbar_inc,
    hist_azimuthal_angle_separation_ttbar_inc, hist_rapidity_difference_ttbar_inc, hist_met_ttbar_inc,
    hist_mT_W_ttbar_inc, hist_m_WW_ttbar_inc, hist_pT_WW_ttbar_inc, hist_leading_lepton_eta_ttbar_inc,
    hist_subleading_lepton_eta_ttbar_inc
)



(hist_leading_lepton_single_top, hist_subleading_lepton_single_top, hist_m_ll_single_top,
 hist_pt_ll_single_top, hist_rapidity_ll_single_top, hist_lepton_separation_single_top,
 hist_azimuthal_angle_separation_single_top, hist_rapidity_difference_single_top, hist_met_single_top,
 hist_mT_W_single_top, hist_m_WW_single_top, hist_pT_WW_single_top, hist_leading_lepton_eta_single_top,
 hist_subleading_lepton_eta_single_top, single_top_efficiency_pre, single_top_efficiency_final) = process_file(
    single_top_background_file_path, hist_leading_lepton_single_top, hist_subleading_lepton_single_top, hist_m_ll_single_top,
    hist_pt_ll_single_top, hist_rapidity_ll_single_top, hist_lepton_separation_single_top,
    hist_azimuthal_angle_separation_single_top, hist_rapidity_difference_single_top, hist_met_single_top,
    hist_mT_W_single_top, hist_m_WW_single_top, hist_pT_WW_single_top, hist_leading_lepton_eta_single_top,
    hist_subleading_lepton_eta_single_top
)




(hist_leading_lepton_w, hist_subleading_lepton_w, hist_m_ll_w,
 hist_pt_ll_w, hist_rapidity_ll_w, hist_lepton_separation_w,
 hist_azimuthal_angle_separation_w, hist_rapidity_difference_w, hist_met_w,
 hist_mT_W_w, hist_m_WW_w, hist_pT_WW_w, hist_leading_lepton_eta_w,
 hist_subleading_lepton_eta_w, w_efficiency_pre, w_efficiency_final) = process_file(
    w_production_background_file_path, hist_leading_lepton_w, hist_subleading_lepton_w, hist_m_ll_w,
    hist_pt_ll_w, hist_rapidity_ll_w, hist_lepton_separation_w,
    hist_azimuthal_angle_separation_w, hist_rapidity_difference_w, hist_met_w,
    hist_mT_W_w, hist_m_WW_w, hist_pT_WW_w, hist_leading_lepton_eta_w,
    hist_subleading_lepton_eta_w
)



(hist_leading_lepton_z, hist_subleading_lepton_z, hist_m_ll_z,
 hist_pt_ll_z, hist_rapidity_ll_z, hist_lepton_separation_z,
 hist_azimuthal_angle_separation_z, hist_rapidity_difference_z, hist_met_z,
 hist_mT_W_z, hist_m_WW_z, hist_pT_WW_z, hist_leading_lepton_eta_z,
 hist_subleading_lepton_eta_z, z_efficiency_pre, z_efficiency_final) = process_file(
    z_production_background_file_path, hist_leading_lepton_z, hist_subleading_lepton_z, hist_m_ll_z,
    hist_pt_ll_z, hist_rapidity_ll_z, hist_lepton_separation_z,
    hist_azimuthal_angle_separation_z, hist_rapidity_difference_z, hist_met_z,
    hist_mT_W_z, hist_m_WW_z, hist_pT_WW_z, hist_leading_lepton_eta_z,
    hist_subleading_lepton_eta_z
)







# Print selection efficiencies for all signals
print("\n=== Signal Selection Efficiencies ===")
for signal_name, efficiencies in signal_efficiencies.items():
    print(f"{signal_name} Selection Efficiency (Pre)  : {efficiencies['efficiency_pre']:.2%}")
    print(f"{signal_name} Selection Efficiency (Final): {efficiencies['efficiency_final']:.2%}")
    print("-" * 50)



# Print background efficiencies
print("\n=== Background Selection Efficiency ===")
print(f"Background Selection Efficiency (Pre)  : {background_efficiency_pre:.2%}")
print(f"Background Selection Efficiency (Final): {background_efficiency_final:.2%}")
print("=" * 50)



print("\n=== Background Selection Efficiencies ===")
print(f"SM WW efficiency (Final)          = {background_efficiency_final:.2%}")
print(f"γγ → τ⁺τ⁻ efficiency (Final)       = {tautau_efficiency_final:.2%}")
print(f"γγ → μ⁺μ⁻ efficiency (Final)       = {mumu_efficiency_final:.2%}")
print(f"γγ → t̄t efficiency (Final)        = {ttbar_efficiency_final:.2%}")
print(f"Inclusive t̄t efficiency (Final)   = {ttbar_inc_efficiency_final:.2%}")
print(f"Single Top efficiency (Final)     = {single_top_efficiency_final:.2%}")
print(f"W Production efficiency (Final)   = {w_efficiency_final:.2%}")
print(f"Z Production efficiency (Final)   = {z_efficiency_final:.2%}")




# Dictionary to store differential cross-sections for each signal
signal_dsigma = {}

# Calculate differential cross-sections for all signals
for signal_name, histograms in signal_histograms.items():
    print(f"Calculating differential cross-sections for {signal_name}...")

    cross_section = signal_cross_sections[signal_name]  # Get cross-section for this signal

    signal_dsigma[signal_name] = {
        # ✅ Leading & Subleading Lepton Kinematics
        "pt_bins_leading_lepton": calculate_dsigma(histograms["hist_leading_lepton"], cross_section, bin_width_pt_leading_lepton),
        "pt_bins_subleading_lepton": calculate_dsigma(histograms["hist_subleading_lepton"], cross_section, bin_width_pt_subleading_lepton),
        "eta_bins_leading_lepton": calculate_dsigma(histograms["hist_leading_lepton_eta"], cross_section, bin_width_eta_leading_lepton),
        "eta_bins_subleading_lepton": calculate_dsigma(histograms["hist_subleading_lepton_eta"], cross_section, bin_width_eta_subleading_lepton),

        # ✅ Dilepton System Kinematics
        "m_ll_bins": calculate_dsigma(histograms["hist_m_ll"], cross_section, bin_width_m_ll),
        "pt_ll_bins": calculate_dsigma(histograms["hist_pt_ll"], cross_section, bin_width_pt_ll),
        "rapidity_ll_bins": calculate_dsigma(histograms["hist_rapidity_ll"], cross_section, bin_width_rapidity_ll),

        # ✅ Angular Separations
        "delta_r_bins": calculate_dsigma(histograms["hist_lepton_separation"], cross_section, bin_width_delta_r),
        "delta_phi_bins": calculate_dsigma(histograms["hist_azimuthal_angle_separation"], cross_section, bin_width_delta_phi),
        "delta_y_bins": calculate_dsigma(histograms["hist_rapidity_difference"], cross_section, bin_width_delta_y),

        # ✅ Missing Transverse Energy
        "met_bins": calculate_dsigma(histograms["hist_met"], cross_section, bin_width_met),

        # ✅ W Boson & Diboson System
        "mT_W_bins": calculate_dsigma(histograms["hist_mT_W"], cross_section, bin_width_mT_W),
        "m_WW_bins": calculate_dsigma(histograms["hist_m_WW"], cross_section, bin_width_m_WW),
        "pT_WW_bins": calculate_dsigma(histograms["hist_pT_WW"], cross_section, bin_width_pT_WW),
    }


# Calculate differential cross-sections for SM WW background
background_dsigma = {
    # ✅ Leading & Subleading Lepton Kinematics
    "pt_bins_leading_lepton": calculate_dsigma(hist_leading_lepton_background, aa_ww_background_cross_section, bin_width_pt_leading_lepton),
    "pt_bins_subleading_lepton": calculate_dsigma(hist_subleading_lepton_background, aa_ww_background_cross_section, bin_width_pt_subleading_lepton),
    "eta_bins_leading_lepton": calculate_dsigma(hist_leading_lepton_eta_background, aa_ww_background_cross_section, bin_width_eta_leading_lepton),
    "eta_bins_subleading_lepton": calculate_dsigma(hist_subleading_lepton_eta_background, aa_ww_background_cross_section, bin_width_eta_subleading_lepton),

    # ✅ Dilepton System Kinematics
    "m_ll_bins": calculate_dsigma(hist_m_ll_background, aa_ww_background_cross_section, bin_width_m_ll),
    "pt_ll_bins": calculate_dsigma(hist_pt_ll_background, aa_ww_background_cross_section, bin_width_pt_ll),
    "rapidity_ll_bins": calculate_dsigma(hist_rapidity_ll_background, aa_ww_background_cross_section, bin_width_rapidity_ll),

    # ✅ Angular Separations
    "delta_r_bins": calculate_dsigma(hist_lepton_separation_background, aa_ww_background_cross_section, bin_width_delta_r),
    "delta_phi_bins": calculate_dsigma(hist_azimuthal_angle_separation_background, aa_ww_background_cross_section, bin_width_delta_phi),
    "delta_y_bins": calculate_dsigma(hist_rapidity_difference_background, aa_ww_background_cross_section, bin_width_delta_y),

    # ✅ Missing Transverse Energy
    "met_bins": calculate_dsigma(hist_met_background, aa_ww_background_cross_section, bin_width_met),

    # ✅ W Boson & Diboson System
    "mT_W_bins": calculate_dsigma(hist_mT_W_background, aa_ww_background_cross_section, bin_width_mT_W),
    "m_WW_bins": calculate_dsigma(hist_m_WW_background, aa_ww_background_cross_section, bin_width_m_WW),
    "pT_WW_bins": calculate_dsigma(hist_pT_WW_background, aa_ww_background_cross_section, bin_width_pT_WW),
}




# Calculate differential cross-sections for ttbar background
ttbar_dsigma = {
    "pt_bins_leading_lepton": calculate_dsigma(hist_leading_lepton_ttbar, aa_ttbar_background_cross_section, bin_width_pt_leading_lepton),
    "pt_bins_subleading_lepton": calculate_dsigma(hist_subleading_lepton_ttbar, aa_ttbar_background_cross_section, bin_width_pt_subleading_lepton),
    "eta_bins_leading_lepton": calculate_dsigma(hist_leading_lepton_eta_ttbar, aa_ttbar_background_cross_section, bin_width_eta_leading_lepton),
    "eta_bins_subleading_lepton": calculate_dsigma(hist_subleading_lepton_eta_ttbar, aa_ttbar_background_cross_section, bin_width_eta_subleading_lepton),
    "m_ll_bins": calculate_dsigma(hist_m_ll_ttbar, aa_ttbar_background_cross_section, bin_width_m_ll),
    "pt_ll_bins": calculate_dsigma(hist_pt_ll_ttbar, aa_ttbar_background_cross_section, bin_width_pt_ll),
    "rapidity_ll_bins": calculate_dsigma(hist_rapidity_ll_ttbar, aa_ttbar_background_cross_section, bin_width_rapidity_ll),
    "delta_r_bins": calculate_dsigma(hist_lepton_separation_ttbar, aa_ttbar_background_cross_section, bin_width_delta_r),
    "delta_phi_bins": calculate_dsigma(hist_azimuthal_angle_separation_ttbar, aa_ttbar_background_cross_section, bin_width_delta_phi),
    "delta_y_bins": calculate_dsigma(hist_rapidity_difference_ttbar, aa_ttbar_background_cross_section, bin_width_delta_y),
    "met_bins": calculate_dsigma(hist_met_ttbar, aa_ttbar_background_cross_section, bin_width_met),
    "mT_W_bins": calculate_dsigma(hist_mT_W_ttbar, aa_ttbar_background_cross_section, bin_width_mT_W),
    "m_WW_bins": calculate_dsigma(hist_m_WW_ttbar, aa_ttbar_background_cross_section, bin_width_m_WW),
    "pT_WW_bins": calculate_dsigma(hist_pT_WW_ttbar, aa_ttbar_background_cross_section, bin_width_pT_WW),
}




# Calculate differential cross-sections for tautau background
tautau_dsigma = {
    "pt_bins_leading_lepton": calculate_dsigma(hist_leading_lepton_tautau, aa_tautau_background_cross_section, bin_width_pt_leading_lepton),
    "pt_bins_subleading_lepton": calculate_dsigma(hist_subleading_lepton_tautau, aa_tautau_background_cross_section, bin_width_pt_subleading_lepton),
    "eta_bins_leading_lepton": calculate_dsigma(hist_leading_lepton_eta_tautau, aa_tautau_background_cross_section, bin_width_eta_leading_lepton),
    "eta_bins_subleading_lepton": calculate_dsigma(hist_subleading_lepton_eta_tautau, aa_tautau_background_cross_section, bin_width_eta_subleading_lepton),
    "m_ll_bins": calculate_dsigma(hist_m_ll_tautau, aa_tautau_background_cross_section, bin_width_m_ll),
    "pt_ll_bins": calculate_dsigma(hist_pt_ll_tautau, aa_tautau_background_cross_section, bin_width_pt_ll),
    "rapidity_ll_bins": calculate_dsigma(hist_rapidity_ll_tautau, aa_tautau_background_cross_section, bin_width_rapidity_ll),
    "delta_r_bins": calculate_dsigma(hist_lepton_separation_tautau, aa_tautau_background_cross_section, bin_width_delta_r),
    "delta_phi_bins": calculate_dsigma(hist_azimuthal_angle_separation_tautau, aa_tautau_background_cross_section, bin_width_delta_phi),
    "delta_y_bins": calculate_dsigma(hist_rapidity_difference_tautau, aa_tautau_background_cross_section, bin_width_delta_y),
    "met_bins": calculate_dsigma(hist_met_tautau, aa_tautau_background_cross_section, bin_width_met),
    "mT_W_bins": calculate_dsigma(hist_mT_W_tautau, aa_tautau_background_cross_section, bin_width_mT_W),
    "m_WW_bins": calculate_dsigma(hist_m_WW_tautau, aa_tautau_background_cross_section, bin_width_m_WW),
    "pT_WW_bins": calculate_dsigma(hist_pT_WW_tautau, aa_tautau_background_cross_section, bin_width_pT_WW),
}





# Calculate differential cross-sections for mumu background
mumu_dsigma = {
    "pt_bins_leading_lepton": calculate_dsigma(hist_leading_lepton_mumu, aa_mumu_background_cross_section, bin_width_pt_leading_lepton),
    "pt_bins_subleading_lepton": calculate_dsigma(hist_subleading_lepton_mumu, aa_mumu_background_cross_section, bin_width_pt_subleading_lepton),
    "eta_bins_leading_lepton": calculate_dsigma(hist_leading_lepton_eta_mumu, aa_mumu_background_cross_section, bin_width_eta_leading_lepton),
    "eta_bins_subleading_lepton": calculate_dsigma(hist_subleading_lepton_eta_mumu, aa_mumu_background_cross_section, bin_width_eta_subleading_lepton),
    "m_ll_bins": calculate_dsigma(hist_m_ll_mumu, aa_mumu_background_cross_section, bin_width_m_ll),
    "pt_ll_bins": calculate_dsigma(hist_pt_ll_mumu, aa_mumu_background_cross_section, bin_width_pt_ll),
    "rapidity_ll_bins": calculate_dsigma(hist_rapidity_ll_mumu, aa_mumu_background_cross_section, bin_width_rapidity_ll),
    "delta_r_bins": calculate_dsigma(hist_lepton_separation_mumu, aa_mumu_background_cross_section, bin_width_delta_r),
    "delta_phi_bins": calculate_dsigma(hist_azimuthal_angle_separation_mumu, aa_mumu_background_cross_section, bin_width_delta_phi),
    "delta_y_bins": calculate_dsigma(hist_rapidity_difference_mumu, aa_mumu_background_cross_section, bin_width_delta_y),
    "met_bins": calculate_dsigma(hist_met_mumu, aa_mumu_background_cross_section, bin_width_met),
    "mT_W_bins": calculate_dsigma(hist_mT_W_mumu, aa_mumu_background_cross_section, bin_width_mT_W),
    "m_WW_bins": calculate_dsigma(hist_m_WW_mumu, aa_mumu_background_cross_section, bin_width_m_WW),
    "pT_WW_bins": calculate_dsigma(hist_pT_WW_mumu, aa_mumu_background_cross_section, bin_width_pT_WW),
}





ttbar_inc_dsigma = {
    "pt_bins_leading_lepton": calculate_dsigma(hist_leading_lepton_ttbar_inc, inclusive_ttbar_background_cross_section, bin_width_pt_leading_lepton),
    "pt_bins_subleading_lepton": calculate_dsigma(hist_subleading_lepton_ttbar_inc, inclusive_ttbar_background_cross_section, bin_width_pt_subleading_lepton),
    "eta_bins_leading_lepton": calculate_dsigma(hist_leading_lepton_eta_ttbar_inc, inclusive_ttbar_background_cross_section, bin_width_eta_leading_lepton),
    "eta_bins_subleading_lepton": calculate_dsigma(hist_subleading_lepton_eta_ttbar_inc, inclusive_ttbar_background_cross_section, bin_width_eta_subleading_lepton),
    "m_ll_bins": calculate_dsigma(hist_m_ll_ttbar_inc, inclusive_ttbar_background_cross_section, bin_width_m_ll),
    "pt_ll_bins": calculate_dsigma(hist_pt_ll_ttbar_inc, inclusive_ttbar_background_cross_section, bin_width_pt_ll),
    "rapidity_ll_bins": calculate_dsigma(hist_rapidity_ll_ttbar_inc, inclusive_ttbar_background_cross_section, bin_width_rapidity_ll),
    "delta_r_bins": calculate_dsigma(hist_lepton_separation_ttbar_inc, inclusive_ttbar_background_cross_section, bin_width_delta_r),
    "delta_phi_bins": calculate_dsigma(hist_azimuthal_angle_separation_ttbar_inc, inclusive_ttbar_background_cross_section, bin_width_delta_phi),
    "delta_y_bins": calculate_dsigma(hist_rapidity_difference_ttbar_inc, inclusive_ttbar_background_cross_section, bin_width_delta_y),
    "met_bins": calculate_dsigma(hist_met_ttbar_inc, inclusive_ttbar_background_cross_section, bin_width_met),
    "mT_W_bins": calculate_dsigma(hist_mT_W_ttbar_inc, inclusive_ttbar_background_cross_section, bin_width_mT_W),
    "m_WW_bins": calculate_dsigma(hist_m_WW_ttbar_inc, inclusive_ttbar_background_cross_section, bin_width_m_WW),
    "pT_WW_bins": calculate_dsigma(hist_pT_WW_ttbar_inc, inclusive_ttbar_background_cross_section, bin_width_pT_WW),
}




single_top_dsigma = {
    "pt_bins_leading_lepton": calculate_dsigma(hist_leading_lepton_single_top, single_top_background_cross_section, bin_width_pt_leading_lepton),
    "pt_bins_subleading_lepton": calculate_dsigma(hist_subleading_lepton_single_top, single_top_background_cross_section, bin_width_pt_subleading_lepton),
    "eta_bins_leading_lepton": calculate_dsigma(hist_leading_lepton_eta_single_top, single_top_background_cross_section, bin_width_eta_leading_lepton),
    "eta_bins_subleading_lepton": calculate_dsigma(hist_subleading_lepton_eta_single_top, single_top_background_cross_section, bin_width_eta_subleading_lepton),
    "m_ll_bins": calculate_dsigma(hist_m_ll_single_top, single_top_background_cross_section, bin_width_m_ll),
    "pt_ll_bins": calculate_dsigma(hist_pt_ll_single_top, single_top_background_cross_section, bin_width_pt_ll),
    "rapidity_ll_bins": calculate_dsigma(hist_rapidity_ll_single_top, single_top_background_cross_section, bin_width_rapidity_ll),
    "delta_r_bins": calculate_dsigma(hist_lepton_separation_single_top, single_top_background_cross_section, bin_width_delta_r),
    "delta_phi_bins": calculate_dsigma(hist_azimuthal_angle_separation_single_top, single_top_background_cross_section, bin_width_delta_phi),
    "delta_y_bins": calculate_dsigma(hist_rapidity_difference_single_top, single_top_background_cross_section, bin_width_delta_y),
    "met_bins": calculate_dsigma(hist_met_single_top, single_top_background_cross_section, bin_width_met),
    "mT_W_bins": calculate_dsigma(hist_mT_W_single_top, single_top_background_cross_section, bin_width_mT_W),
    "m_WW_bins": calculate_dsigma(hist_m_WW_single_top, single_top_background_cross_section, bin_width_m_WW),
    "pT_WW_bins": calculate_dsigma(hist_pT_WW_single_top, single_top_background_cross_section, bin_width_pT_WW),
}




w_dsigma = {
    "pt_bins_leading_lepton": calculate_dsigma(hist_leading_lepton_w, w_production_background_cross_section, bin_width_pt_leading_lepton),
    "pt_bins_subleading_lepton": calculate_dsigma(hist_subleading_lepton_w, w_production_background_cross_section, bin_width_pt_subleading_lepton),
    "eta_bins_leading_lepton": calculate_dsigma(hist_leading_lepton_eta_w, w_production_background_cross_section, bin_width_eta_leading_lepton),
    "eta_bins_subleading_lepton": calculate_dsigma(hist_subleading_lepton_eta_w, w_production_background_cross_section, bin_width_eta_subleading_lepton),
    "m_ll_bins": calculate_dsigma(hist_m_ll_w, w_production_background_cross_section, bin_width_m_ll),
    "pt_ll_bins": calculate_dsigma(hist_pt_ll_w, w_production_background_cross_section, bin_width_pt_ll),
    "rapidity_ll_bins": calculate_dsigma(hist_rapidity_ll_w, w_production_background_cross_section, bin_width_rapidity_ll),
    "delta_r_bins": calculate_dsigma(hist_lepton_separation_w, w_production_background_cross_section, bin_width_delta_r),
    "delta_phi_bins": calculate_dsigma(hist_azimuthal_angle_separation_w, w_production_background_cross_section, bin_width_delta_phi),
    "delta_y_bins": calculate_dsigma(hist_rapidity_difference_w, w_production_background_cross_section, bin_width_delta_y),
    "met_bins": calculate_dsigma(hist_met_w, w_production_background_cross_section, bin_width_met),
    "mT_W_bins": calculate_dsigma(hist_mT_W_w, w_production_background_cross_section, bin_width_mT_W),
    "m_WW_bins": calculate_dsigma(hist_m_WW_w, w_production_background_cross_section, bin_width_m_WW),
    "pT_WW_bins": calculate_dsigma(hist_pT_WW_w, w_production_background_cross_section, bin_width_pT_WW),
}





z_dsigma = {
    "pt_bins_leading_lepton": calculate_dsigma(hist_leading_lepton_z, z_production_background_cross_section, bin_width_pt_leading_lepton),
    "pt_bins_subleading_lepton": calculate_dsigma(hist_subleading_lepton_z, z_production_background_cross_section, bin_width_pt_subleading_lepton),
    "eta_bins_leading_lepton": calculate_dsigma(hist_leading_lepton_eta_z, z_production_background_cross_section, bin_width_eta_leading_lepton),
    "eta_bins_subleading_lepton": calculate_dsigma(hist_subleading_lepton_eta_z, z_production_background_cross_section, bin_width_eta_subleading_lepton),
    "m_ll_bins": calculate_dsigma(hist_m_ll_z, z_production_background_cross_section, bin_width_m_ll),
    "pt_ll_bins": calculate_dsigma(hist_pt_ll_z, z_production_background_cross_section, bin_width_pt_ll),
    "rapidity_ll_bins": calculate_dsigma(hist_rapidity_ll_z, z_production_background_cross_section, bin_width_rapidity_ll),
    "delta_r_bins": calculate_dsigma(hist_lepton_separation_z, z_production_background_cross_section, bin_width_delta_r),
    "delta_phi_bins": calculate_dsigma(hist_azimuthal_angle_separation_z, z_production_background_cross_section, bin_width_delta_phi),
    "delta_y_bins": calculate_dsigma(hist_rapidity_difference_z, z_production_background_cross_section, bin_width_delta_y),
    "met_bins": calculate_dsigma(hist_met_z, z_production_background_cross_section, bin_width_met),
    "mT_W_bins": calculate_dsigma(hist_mT_W_z, z_production_background_cross_section, bin_width_mT_W),
    "m_WW_bins": calculate_dsigma(hist_m_WW_z, z_production_background_cross_section, bin_width_m_WW),
    "pT_WW_bins": calculate_dsigma(hist_pT_WW_z, z_production_background_cross_section, bin_width_pT_WW),
}






#=========================================================================
#=========================================================================


plt.figure(figsize=(11, 12))  # Create a new figure for the leading lepton pT plot
plt.subplots_adjust(left=0.15, right=0.95, bottom=0.12, top=0.95)

# Define colors for each signal
signal_colors = {
    "$FM_{0} / \Lambda^4$": "red",
    "$FM_{1} / \Lambda^4$": "purple",
    "$FM_{2} / \Lambda^4$": "green",
    "$FM_{3} / \Lambda^4$": "orange"
}





# ===================================================
# ✅ Leading Lepton \( p_T \) Differential Cross-Section
# ===================================================

for signal_name, dsigma_data in signal_dsigma.items():
    pt_bins, dsigma = dsigma_data["pt_bins_leading_lepton"]
    plt.step(pt_bins, dsigma, where="mid", alpha=0.7,
             label=f"LHeC@1.2 TeV : Signal ($W^+ W^-$) [{signal_name}]",
             color=signal_colors.get(signal_name, "black"), linewidth=3)

# ✅ Plot Background for Leading Lepton
pt_bins_background, dsigma_background = background_dsigma["pt_bins_leading_lepton"]
plt.step(pt_bins_background, dsigma_background, where="mid", alpha=0.7,
         label=r"LHeC@1.2 TeV : SM ($\gamma\gamma \to W^+ W^-$)", color="blue", linewidth=3)

pt_bins_ttbar, dsigma_ttbar = ttbar_dsigma["pt_bins_leading_lepton"]
plt.step(pt_bins_ttbar, dsigma_ttbar, where="mid", alpha=0.7, label=r"LHeC@1.2 TeV : $\gamma\gamma \to t\bar{t} (\times 10^{3})$"
, color="purple", linewidth=3)


# ✅ tautau Background
pt_bins_tautau, dsigma_tautau = tautau_dsigma["pt_bins_leading_lepton"]
plt.step(pt_bins_tautau, dsigma_tautau, where="mid", alpha=0.7,
         label=r"LHeC@1.2 TeV : $\gamma\gamma \to \tau^+\tau^-$", color="orange", linewidth=3)

# ✅ mumu Background
pt_bins_mumu, dsigma_mumu = mumu_dsigma["pt_bins_leading_lepton"]
plt.step(pt_bins_mumu, dsigma_mumu, where="mid", alpha=0.7,
         label=r"LHeC@1.2 TeV : $\gamma\gamma \to \mu^+\mu^-$", color="brown", linewidth=3)


# ✅ Inclusive tt̄ Background
pt_bins_ttbar_inc, dsigma_ttbar_inc = ttbar_inc_dsigma["pt_bins_leading_lepton"]
plt.step(pt_bins_ttbar_inc, dsigma_ttbar_inc, where="mid", alpha=0.7,
         label=r"LHeC@1.2 TeV : Inclusive $t\bar{t}$", color="orchid", linewidth=3)

# ✅ Single Top Background
pt_bins_single_top, dsigma_single_top = single_top_dsigma["pt_bins_leading_lepton"]
plt.step(pt_bins_single_top, dsigma_single_top, where="mid", alpha=0.7,
         label=r"LHeC@1.2 TeV : Single Top", color="teal", linewidth=3)

# ✅ W Production Background
pt_bins_w, dsigma_w = w_dsigma["pt_bins_leading_lepton"]
plt.step(pt_bins_w, dsigma_w, where="mid", alpha=0.7,
         label=r"LHeC@1.2 TeV : W Production", color="darkgreen", linewidth=3)

# ✅ Z Production Background
pt_bins_z, dsigma_z = z_dsigma["pt_bins_leading_lepton"]
plt.step(pt_bins_z, dsigma_z, where="mid", alpha=0.7,
         label=r"LHeC@1.2 TeV : Z Production", color="darkred", linewidth=3)



# Set labels and title
plt.xlabel(r"$p_T^{\mathrm{leading}~\ell} \ \mathrm{[GeV]}$")
plt.ylabel(r"$\frac{d\sigma}{dp_T^{\mathrm{leading}~\ell}} \ \mathrm{[pb/GeV]}$")
plt.title(r"Delphes simulation : $e^- p \to e^- W^+ W^- p \to e^- \ell^+ \nu_{\ell} \ell^- \bar{\nu}_{\ell} p$ : LHeC@1.2 TeV", fontsize=20)
plt.yscale("log")

# Grid, legend, and layout adjustments
plt.legend()
plt.grid(True, linestyle="--", alpha=0.6)
plt.tight_layout()
plt.ylim(0.00001, 0.001)

# Save the plot
plt.savefig("/home/hamzeh-khanpour/Documents/GitHub/LHeC_Fast_Simulation/Pythia8_Delphes_fully_leptonic_allsignal_bkgs/differential_cross_section_leading_lepton_pt_allFMsignal_allbkgs.pdf", dpi=600)

# Show the plot
plt.show()







# ===================================================
# ✅ Subleading Lepton \( p_T \) Differential Cross-Section
# ===================================================

plt.figure(figsize=(11, 12))  # Create a new figure for the subleading lepton pT plot
plt.subplots_adjust(left=0.15, right=0.95, bottom=0.12, top=0.95)

# Loop through all signals and plot their differential cross-section for subleading lepton
for signal_name, dsigma_data in signal_dsigma.items():
    pt_bins, dsigma = dsigma_data["pt_bins_subleading_lepton"]
    plt.step(pt_bins, dsigma, where="mid", alpha=0.7,
             label=f"LHeC@1.2 TeV : Signal ($W^+ W^-$) [{signal_name}]",
             color=signal_colors.get(signal_name, "black"), linewidth=3)

# ✅ Plot Background for Subleading Lepton
pt_bins_background, dsigma_background = background_dsigma["pt_bins_subleading_lepton"]
plt.step(pt_bins_background, dsigma_background, where="mid", alpha=0.7,
         label=r"LHeC@1.2 TeV : SM ($\gamma\gamma \to W^+ W^-$)", color="blue", linewidth=3)

# ✅ Plot ttbar Background
pt_bins_ttbar, dsigma_ttbar = ttbar_dsigma["pt_bins_subleading_lepton"]
plt.step(pt_bins_ttbar, dsigma_ttbar, where="mid", alpha=0.7,
         label=r"LHeC@1.2 TeV : $\gamma\gamma \to t\bar{t} (\times 10^{3})$", color="purple", linewidth=3)


# ✅ tautau Background
pt_bins_tautau, dsigma_tautau = tautau_dsigma["pt_bins_subleading_lepton"]
plt.step(pt_bins_tautau, dsigma_tautau, where="mid", alpha=0.7,
         label=r"LHeC@1.2 TeV : $\gamma\gamma \to \tau^+\tau^-$", color="orange", linewidth=3)

# ✅ mumu Background
pt_bins_mumu, dsigma_mumu = mumu_dsigma["pt_bins_subleading_lepton"]
plt.step(pt_bins_mumu, dsigma_mumu, where="mid", alpha=0.7,
         label=r"LHeC@1.2 TeV : $\gamma\gamma \to \mu^+\mu^-$", color="brown", linewidth=3)



# Set labels and title
plt.xlabel(r"$p_T^{\mathrm{subleading}~\ell} \ \mathrm{[GeV]}$")
plt.ylabel(r"$\frac{d\sigma}{dp_T^{\mathrm{subleading}~\ell}} \ \mathrm{[pb/GeV]}$")
plt.title(r"Delphes simulation : $e^- p \to e^- W^+ W^- p \to e^- \ell^+ \nu_{\ell} \ell^- \bar{\nu}_{\ell} p$ : LHeC@1.2 TeV", fontsize=20)
plt.yscale("log")

# Grid, legend, and layout adjustments
plt.legend()
plt.grid(True, linestyle="--", alpha=0.6)
plt.tight_layout()
plt.ylim(0.00001, 0.001)

# Save the plot
plt.savefig("/home/hamzeh-khanpour/Documents/GitHub/LHeC_Fast_Simulation/Pythia8_Delphes_fully_leptonic_allsignal_bkgs/differential_cross_section_subleading_lepton_pt_allFMsignal_allbkgs.pdf", dpi=600)

# Show the plot
plt.show()








# ===================================================
# ✅ Dilepton Invariant Mass \( M_{\ell\ell} \) Differential Cross-Section
# ===================================================

plt.figure(figsize=(11, 12))  # Create a new figure for M_ll plot
plt.subplots_adjust(left=0.15, right=0.95, bottom=0.12, top=0.95)

# Loop through all signals and plot their differential cross-section for M_ll
for signal_name, dsigma_data in signal_dsigma.items():
    m_ll_bins, dsigma = dsigma_data["m_ll_bins"]
    plt.step(m_ll_bins, dsigma, where="mid", alpha=0.7,
             label=f"LHeC@1.2 TeV : Signal ($W^+ W^-$) [{signal_name}]",
             color=signal_colors.get(signal_name, "black"), linewidth=3)

# ✅ Plot Background for Dilepton Invariant Mass
m_ll_bins_background, dsigma_background = background_dsigma["m_ll_bins"]
plt.step(m_ll_bins_background, dsigma_background, where="mid", alpha=0.7,
         label=r"LHeC@1.2 TeV : SM ($\gamma\gamma \to W^+ W^-$)", color="blue", linewidth=3)


# ✅ Plot ttbar Background
m_ll_bins_ttbar, dsigma_ttbar = ttbar_dsigma["m_ll_bins"]
plt.step(m_ll_bins_ttbar, dsigma_ttbar, where="mid", alpha=0.7,
         label=r"LHeC@1.2 TeV : $\gamma\gamma \to t\bar{t} (\times 10^{3})$", color="purple", linewidth=3)


# ✅ Plot tautau Background
m_ll_bins_tautau, dsigma_tautau = tautau_dsigma["m_ll_bins"]
plt.step(m_ll_bins_tautau, dsigma_tautau, where="mid", alpha=0.7,
         label=r"LHeC@1.2 TeV : $\gamma\gamma \to \tau^+\tau^-$", color="orange", linewidth=3)

# ✅ Plot mumu Background
m_ll_bins_mumu, dsigma_mumu = mumu_dsigma["m_ll_bins"]
plt.step(m_ll_bins_mumu, dsigma_mumu, where="mid", alpha=0.7,
         label=r"LHeC@1.2 TeV : $\gamma\gamma \to \mu^+\mu^-$", color="brown", linewidth=3)

# ✅ Inclusive tt̄ Background
m_ll_bins_ttbar_inc, dsigma_ttbar_inc = ttbar_inc_dsigma["m_ll_bins"]
plt.step(m_ll_bins_ttbar_inc, dsigma_ttbar_inc, where="mid", alpha=0.7,
         label=r"LHeC@1.2 TeV : Inclusive $t\bar{t}$", color="orchid", linewidth=3)

# ✅ Single Top Background
m_ll_bins_single_top, dsigma_single_top = single_top_dsigma["m_ll_bins"]
plt.step(m_ll_bins_single_top, dsigma_single_top, where="mid", alpha=0.7,
         label=r"LHeC@1.2 TeV : Single Top", color="teal", linewidth=3)

# ✅ W Production Background
m_ll_bins_w, dsigma_w = w_dsigma["m_ll_bins"]
plt.step(m_ll_bins_w, dsigma_w, where="mid", alpha=0.7,
         label=r"LHeC@1.2 TeV : W Production", color="darkgreen", linewidth=3)

# ✅ Z Production Background
m_ll_bins_z, dsigma_z = z_dsigma["m_ll_bins"]
plt.step(m_ll_bins_z, dsigma_z, where="mid", alpha=0.7,
         label=r"LHeC@1.2 TeV : Z Production", color="darkred", linewidth=3)


# Set labels and title
plt.xlabel(r"$M_{\ell\ell} \ \mathrm{[GeV]}$")
plt.ylabel(r"$\frac{d\sigma}{dM_{\ell\ell}} \ \mathrm{[pb/GeV]}$")
plt.title(r"Delphes simulation : $e^- p \to e^- W^+ W^- p \to e^- \ell^+ \nu_{\ell} \ell^- \bar{\nu}_{\ell} p$ : LHeC@1.2 TeV", fontsize=20)
plt.yscale("log")

# Grid, legend, and layout adjustments
plt.legend()
plt.grid(True, linestyle="--", alpha=0.6)
plt.tight_layout()
plt.ylim(0.00001, 0.001)

# Save the plot
plt.savefig("/home/hamzeh-khanpour/Documents/GitHub/LHeC_Fast_Simulation/Pythia8_Delphes_fully_leptonic_allsignal_bkgs/differential_cross_section_m_ll_allFMsignal_allbkgs.pdf", dpi=600)

# Show the plot
plt.show()




# ===================================================
# ✅ Dilepton Transverse Momentum \( p_T^{\ell\ell} \) Differential Cross-Section
# ===================================================

plt.figure(figsize=(11, 12))  # Create a new figure for pT_ll plot
plt.subplots_adjust(left=0.15, right=0.95, bottom=0.12, top=0.95)

# Loop through all signals and plot their differential cross-section for pT_ll
for signal_name, dsigma_data in signal_dsigma.items():
    pt_ll_bins, dsigma = dsigma_data["pt_ll_bins"]
    plt.step(pt_ll_bins, dsigma, where="mid", alpha=0.7,
             label=f"LHeC@1.2 TeV : Signal ($W^+ W^-$) [{signal_name}]",
             color=signal_colors.get(signal_name, "black"), linewidth=3)

# ✅ Plot Background for Dilepton Transverse Momentum
pt_ll_bins_background, dsigma_background = background_dsigma["pt_ll_bins"]
plt.step(pt_ll_bins_background, dsigma_background, where="mid", alpha=0.7,
         label=r"LHeC@1.2 TeV : SM ($\gamma\gamma \to W^+ W^-$)", color="blue", linewidth=3)

# ✅ Plot tt̄ Background
pt_ll_bins_ttbar, dsigma_ttbar = ttbar_dsigma["pt_ll_bins"]
plt.step(pt_ll_bins_ttbar, dsigma_ttbar, where="mid", alpha=0.7,
         label=r"LHeC@1.2 TeV : $\gamma\gamma \to t\bar{t} (\times 10^{3})$", color="purple", linewidth=3)


# ✅ Plot τ⁺τ⁻ Background
pt_ll_bins_tautau, dsigma_tautau = tautau_dsigma["pt_ll_bins"]
plt.step(pt_ll_bins_tautau, dsigma_tautau, where="mid", alpha=0.7,
         label=r"LHeC@1.2 TeV : $\gamma\gamma \to \tau^+\tau^-$", color="orange", linewidth=3)

# ✅ Plot μ⁺μ⁻ Background
pt_ll_bins_mumu, dsigma_mumu = mumu_dsigma["pt_ll_bins"]
plt.step(pt_ll_bins_mumu, dsigma_mumu, where="mid", alpha=0.7,
         label=r"LHeC@1.2 TeV : $\gamma\gamma \to \mu^+\mu^-$", color="brown", linewidth=3)

# ✅ Inclusive tt̄ Background
pt_ll_bins_ttbar_inc, dsigma_ttbar_inc = ttbar_inc_dsigma["pt_ll_bins"]
plt.step(pt_ll_bins_ttbar_inc, dsigma_ttbar_inc, where="mid", alpha=0.7,
         label=r"LHeC@1.2 TeV : Inclusive $t\bar{t}$", color="orchid", linewidth=3)

# ✅ Single Top Background
pt_ll_bins_single_top, dsigma_single_top = single_top_dsigma["pt_ll_bins"]
plt.step(pt_ll_bins_single_top, dsigma_single_top, where="mid", alpha=0.7,
         label=r"LHeC@1.2 TeV : Single Top", color="teal", linewidth=3)

# ✅ W Production Background
pt_ll_bins_w, dsigma_w = w_dsigma["pt_ll_bins"]
plt.step(pt_ll_bins_w, dsigma_w, where="mid", alpha=0.7,
         label=r"LHeC@1.2 TeV : W Production", color="darkgreen", linewidth=3)

# ✅ Z Production Background
pt_ll_bins_z, dsigma_z = z_dsigma["pt_ll_bins"]
plt.step(pt_ll_bins_z, dsigma_z, where="mid", alpha=0.7,
         label=r"LHeC@1.2 TeV : Z Production", color="darkred", linewidth=3)


# Set labels and title
plt.xlabel(r"$p_T^{\ell\ell} \ \mathrm{[GeV]}$")
plt.ylabel(r"$\frac{d\sigma}{dp_T^{\ell\ell}} \ \mathrm{[pb/GeV]}$")
plt.title(r"Delphes simulation : $e^- p \to e^- W^+ W^- p \to e^- \ell^+ \nu_{\ell} \ell^- \bar{\nu}_{\ell} p$ : LHeC@1.2 TeV", fontsize=20)
plt.yscale("log")

# Grid, legend, and layout adjustments
plt.legend()
plt.grid(True, linestyle="--", alpha=0.6)
plt.tight_layout()
plt.ylim(0.00001, 0.001)

# Save the plot
plt.savefig("/home/hamzeh-khanpour/Documents/GitHub/LHeC_Fast_Simulation/Pythia8_Delphes_fully_leptonic_allsignal_bkgs/differential_cross_section_pt_ll_allFMsignal_allbkgs.pdf", dpi=600)

# Show the plot
plt.show()




# ===================================================
# ✅ Dilepton Rapidity \( Y_{\ell\ell} \) Differential Cross-Section
# ===================================================

plt.figure(figsize=(11, 12))  # Create a new figure for rapidity_ll plot
plt.subplots_adjust(left=0.15, right=0.95, bottom=0.12, top=0.95)

# Loop through all signals and plot their differential cross-section for rapidity_ll
for signal_name, dsigma_data in signal_dsigma.items():
    rapidity_ll_bins, dsigma = dsigma_data["rapidity_ll_bins"]
    plt.step(rapidity_ll_bins, dsigma, where="mid", alpha=0.7,
             label=f"LHeC@1.2 TeV : Signal ($W^+ W^-$) [{signal_name}]",
             color=signal_colors.get(signal_name, "black"), linewidth=3)

# ✅ Plot Background for Dilepton Rapidity
rapidity_ll_bins_background, dsigma_background = background_dsigma["rapidity_ll_bins"]
plt.step(rapidity_ll_bins_background, dsigma_background, where="mid", alpha=0.7,
         label=r"LHeC@1.2 TeV : SM ($\gamma\gamma \to W^+ W^-$)", color="blue", linewidth=3)


# ✅ Plot tt̄ Background
rapidity_ll_bins_ttbar, dsigma_ttbar = ttbar_dsigma["rapidity_ll_bins"]
plt.step(rapidity_ll_bins_ttbar, dsigma_ttbar, where="mid", alpha=0.7,
         label=r"LHeC@1.2 TeV : $\gamma\gamma \to t\bar{t} (\times 10^{3})$", color="purple", linewidth=3)


# ✅ Plot τ⁺τ⁻ Background
rapidity_ll_bins_tautau, dsigma_tautau = tautau_dsigma["rapidity_ll_bins"]
plt.step(rapidity_ll_bins_tautau, dsigma_tautau, where="mid", alpha=0.7,
         label=r"LHeC@1.2 TeV : $\gamma\gamma \to \tau^+\tau^-$", color="orange", linewidth=3)

# ✅ Plot μ⁺μ⁻ Background
rapidity_ll_bins_mumu, dsigma_mumu = mumu_dsigma["rapidity_ll_bins"]
plt.step(rapidity_ll_bins_mumu, dsigma_mumu, where="mid", alpha=0.7,
         label=r"LHeC@1.2 TeV : $\gamma\gamma \to \mu^+\mu^-$", color="brown", linewidth=3)


# ✅ Inclusive tt̄ Background
rapidity_ll_bins_ttbar_inc, dsigma_ttbar_inc = ttbar_inc_dsigma["rapidity_ll_bins"]
plt.step(rapidity_ll_bins_ttbar_inc, dsigma_ttbar_inc, where="mid", alpha=0.7,
         label=r"LHeC@1.2 TeV : Inclusive $t\bar{t}$", color="orchid", linewidth=3)

# ✅ Single Top Background
rapidity_ll_bins_single_top, dsigma_single_top = single_top_dsigma["rapidity_ll_bins"]
plt.step(rapidity_ll_bins_single_top, dsigma_single_top, where="mid", alpha=0.7,
         label=r"LHeC@1.2 TeV : Single Top", color="teal", linewidth=3)

# ✅ W Production Background
rapidity_ll_bins_w, dsigma_w = w_dsigma["rapidity_ll_bins"]
plt.step(rapidity_ll_bins_w, dsigma_w, where="mid", alpha=0.7,
         label=r"LHeC@1.2 TeV : W Production", color="darkgreen", linewidth=3)

# ✅ Z Production Background
rapidity_ll_bins_z, dsigma_z = z_dsigma["rapidity_ll_bins"]
plt.step(rapidity_ll_bins_z, dsigma_z, where="mid", alpha=0.7,
         label=r"LHeC@1.2 TeV : Z Production", color="darkred", linewidth=3)


# Set labels and title
plt.xlabel(r"$Y_{\ell\ell}$")
plt.ylabel(r"$\frac{d\sigma}{dY_{\ell\ell}} \ \mathrm{[pb]}$")
plt.title(r"Delphes simulation : $e^- p \to e^- W^+ W^- p \to e^- \ell^+ \nu_{\ell} \ell^- \bar{\nu}_{\ell} p$ : LHeC@1.2 TeV", fontsize=20)
plt.yscale("log")

# Grid, legend, and layout adjustments
plt.legend()
plt.grid(True, linestyle="--", alpha=0.6)
plt.tight_layout()
plt.ylim(0.00001, 1.0)

# Save the plot
plt.savefig("/home/hamzeh-khanpour/Documents/GitHub/LHeC_Fast_Simulation/Pythia8_Delphes_fully_leptonic_allsignal_bkgs/differential_cross_section_rapidity_ll_allFMsignal_allbkgs.pdf", dpi=600)

# Show the plot
plt.show()




# ===================================================
# ✅ Lepton Separation \( \Delta R(\ell_1, \ell_2) \) Differential Cross-Section
# ===================================================

plt.figure(figsize=(11, 12))  # Create a new figure for Delta R plot
plt.subplots_adjust(left=0.15, right=0.95, bottom=0.12, top=0.95)

# Loop through all signals and plot their differential cross-section for lepton separation
for signal_name, dsigma_data in signal_dsigma.items():
    delta_r_bins, dsigma = dsigma_data["delta_r_bins"]
    plt.step(delta_r_bins, dsigma, where="mid", alpha=0.7,
             label=f"LHeC@1.2 TeV : Signal ($W^+ W^-$) [{signal_name}]",
             color=signal_colors.get(signal_name, "black"), linewidth=3)

# ✅ Plot Background for Lepton Separation
delta_r_bins_background, dsigma_background = background_dsigma["delta_r_bins"]
plt.step(delta_r_bins_background, dsigma_background, where="mid", alpha=0.7,
         label=r"LHeC@1.2 TeV : SM ($\gamma\gamma \to W^+ W^-$)", color="blue", linewidth=3)


# ✅ Plot tt̄ Background
delta_r_bins_ttbar, dsigma_ttbar = ttbar_dsigma["delta_r_bins"]
plt.step(delta_r_bins_ttbar, dsigma_ttbar, where="mid", alpha=0.7,
         label=r"LHeC@1.2 TeV : $\gamma\gamma \to t\bar{t} (\times 10^{3})$", color="purple", linewidth=3)


# ✅ Plot τ⁺τ⁻ Background
delta_r_bins_tautau, dsigma_tautau = tautau_dsigma["delta_r_bins"]
plt.step(delta_r_bins_tautau, dsigma_tautau, where="mid", alpha=0.7,
         label=r"LHeC@1.2 TeV : $\gamma\gamma \to \tau^+\tau^-$", color="orange", linewidth=3)

# ✅ Plot μ⁺μ⁻ Background
delta_r_bins_mumu, dsigma_mumu = mumu_dsigma["delta_r_bins"]
plt.step(delta_r_bins_mumu, dsigma_mumu, where="mid", alpha=0.7,
         label=r"LHeC@1.2 TeV : $\gamma\gamma \to \mu^+\mu^-$", color="brown", linewidth=3)

# ✅ Inclusive tt̄ Background
delta_r_bins_ttbar_inc, dsigma_ttbar_inc = ttbar_inc_dsigma["delta_r_bins"]
plt.step(delta_r_bins_ttbar_inc, dsigma_ttbar_inc, where="mid", alpha=0.7,
         label=r"LHeC@1.2 TeV : Inclusive $t\bar{t}$", color="orchid", linewidth=3)

# ✅ Single Top Background
delta_r_bins_single_top, dsigma_single_top = single_top_dsigma["delta_r_bins"]
plt.step(delta_r_bins_single_top, dsigma_single_top, where="mid", alpha=0.7,
         label=r"LHeC@1.2 TeV : Single Top", color="teal", linewidth=3)

# ✅ W Production Background
delta_r_bins_w, dsigma_w = w_dsigma["delta_r_bins"]
plt.step(delta_r_bins_w, dsigma_w, where="mid", alpha=0.7,
         label=r"LHeC@1.2 TeV : W Production", color="darkgreen", linewidth=3)

# ✅ Z Production Background
delta_r_bins_z, dsigma_z = z_dsigma["delta_r_bins"]
plt.step(delta_r_bins_z, dsigma_z, where="mid", alpha=0.7,
         label=r"LHeC@1.2 TeV : Z Production", color="darkred", linewidth=3)


# Set labels and title
plt.xlabel(r"$\Delta R(\ell_1, \ell_2)$")
plt.ylabel(r"$\frac{d\sigma}{d\Delta R} \ \mathrm{[pb]}$")
plt.title(r"Delphes simulation : $e^- p \to e^- W^+ W^- p \to e^- \ell^+ \nu_{\ell} \ell^- \bar{\nu}_{\ell} p$ : LHeC@1.2 TeV", fontsize=20)
plt.yscale("log")

# Grid, legend, and layout adjustments
plt.legend()
plt.grid(True, linestyle="--", alpha=0.6)
plt.tight_layout()
plt.ylim(0.00001, 1.0)

# Save the plot
plt.savefig("/home/hamzeh-khanpour/Documents/GitHub/LHeC_Fast_Simulation/Pythia8_Delphes_fully_leptonic_allsignal_bkgs/differential_cross_section_deltaR_ll_allFMsignal_allbkgs.pdf", dpi=600)

# Show the plot
plt.show()



# ===================================================
# ✅ Azimuthal Angle Separation \( \Delta \phi(\ell_1, \ell_2) \) Differential Cross-Section
# ===================================================

plt.figure(figsize=(11, 12))  # Create a new figure for Delta Phi plot
plt.subplots_adjust(left=0.15, right=0.95, bottom=0.12, top=0.95)

# Loop through all signals and plot their differential cross-section for azimuthal angle separation
for signal_name, dsigma_data in signal_dsigma.items():
    delta_phi_bins, dsigma = dsigma_data["delta_phi_bins"]
    plt.step(delta_phi_bins, dsigma, where="mid", alpha=0.7,
             label=f"LHeC@1.2 TeV : Signal ($W^+ W^-$) [{signal_name}]",
             color=signal_colors.get(signal_name, "black"), linewidth=3)

# ✅ Plot Background for Azimuthal Angle Separation
delta_phi_bins_background, dsigma_background = background_dsigma["delta_phi_bins"]
plt.step(delta_phi_bins_background, dsigma_background, where="mid", alpha=0.7,
         label=r"LHeC@1.2 TeV : SM ($\gamma\gamma \to W^+ W^-$)", color="blue", linewidth=3)


# ✅ Plot tt̄ Background
delta_phi_bins_ttbar, dsigma_ttbar = ttbar_dsigma["delta_phi_bins"]
plt.step(delta_phi_bins_ttbar, dsigma_ttbar, where="mid", alpha=0.7,
         label=r"LHeC@1.2 TeV : $\gamma\gamma \to t\bar{t} (\times 10^{3})$", color="purple", linewidth=3)

# ✅ Plot τ⁺τ⁻ Background
delta_phi_bins_tautau, dsigma_tautau = tautau_dsigma["delta_phi_bins"]
plt.step(delta_phi_bins_tautau, dsigma_tautau, where="mid", alpha=0.7,
         label=r"LHeC@1.2 TeV : $\gamma\gamma \to \tau^+\tau^-$", color="orange", linewidth=3)

# ✅ Plot μ⁺μ⁻ Background
delta_phi_bins_mumu, dsigma_mumu = mumu_dsigma["delta_phi_bins"]
plt.step(delta_phi_bins_mumu, dsigma_mumu, where="mid", alpha=0.7,
         label=r"LHeC@1.2 TeV : $\gamma\gamma \to \mu^+\mu^-$", color="brown", linewidth=3)

# ✅ Inclusive tt̄ Background
delta_phi_bins_ttbar_inc, dsigma_ttbar_inc = ttbar_inc_dsigma["delta_phi_bins"]
plt.step(delta_phi_bins_ttbar_inc, dsigma_ttbar_inc, where="mid", alpha=0.7,
         label=r"LHeC@1.2 TeV : Inclusive $t\bar{t}$", color="orchid", linewidth=3)

# ✅ Single Top Background
delta_phi_bins_single_top, dsigma_single_top = single_top_dsigma["delta_phi_bins"]
plt.step(delta_phi_bins_single_top, dsigma_single_top, where="mid", alpha=0.7,
         label=r"LHeC@1.2 TeV : Single Top", color="teal", linewidth=3)

# ✅ W Production Background
delta_phi_bins_w, dsigma_w = w_dsigma["delta_phi_bins"]
plt.step(delta_phi_bins_w, dsigma_w, where="mid", alpha=0.7,
         label=r"LHeC@1.2 TeV : W Production", color="darkgreen", linewidth=3)

# ✅ Z Production Background
delta_phi_bins_z, dsigma_z = z_dsigma["delta_phi_bins"]
plt.step(delta_phi_bins_z, dsigma_z, where="mid", alpha=0.7,
         label=r"LHeC@1.2 TeV : Z Production", color="darkred", linewidth=3)


# Set labels and title
plt.xlabel(r"$\Delta \phi(\ell_1, \ell_2)$")
plt.ylabel(r"$\frac{d\sigma}{d\Delta \phi} \ \mathrm{[pb]}$")
plt.title(r"Delphes simulation : $e^- p \to e^- W^+ W^- p \to e^- \ell^+ \nu_{\ell} \ell^- \bar{\nu}_{\ell} p$ : LHeC@1.2 TeV", fontsize=20)
plt.yscale("log")

# Grid, legend, and layout adjustments
plt.legend()
plt.grid(True, linestyle="--", alpha=0.6)
plt.tight_layout()
plt.ylim(0.00001, 10.0)

# Save the plot
plt.savefig("/home/hamzeh-khanpour/Documents/GitHub/LHeC_Fast_Simulation/Pythia8_Delphes_fully_leptonic_allsignal_bkgs/differential_cross_section_deltaPhi_ll_allFMsignal_allbkgs.pdf", dpi=600)

# Show the plot
plt.show()




# ===================================================
# ✅ Rapidity Difference \( \Delta y(\ell_1, \ell_2) \) Differential Cross-Section
# ===================================================

plt.figure(figsize=(11, 12))  # Create a new figure for rapidity difference Δy plot
plt.subplots_adjust(left=0.15, right=0.95, bottom=0.12, top=0.95)

# Loop through all signals and plot their differential cross-section for rapidity difference Δy
for signal_name, dsigma_data in signal_dsigma.items():
    delta_y_bins, dsigma = dsigma_data["delta_y_bins"]
    plt.step(delta_y_bins, dsigma, where="mid", alpha=0.7,
             label=f"LHeC@1.2 TeV : Signal ($W^+ W^-$) [{signal_name}]",
             color=signal_colors.get(signal_name, "black"), linewidth=3)

# ✅ Plot Background for Rapidity Difference
delta_y_bins_background, dsigma_background = background_dsigma["delta_y_bins"]
plt.step(delta_y_bins_background, dsigma_background, where="mid", alpha=0.7,
         label=r"LHeC@1.2 TeV : SM ($\gamma\gamma \to W^+ W^-$)", color="blue", linewidth=3)


# ✅ Plot tt̄ Background
delta_y_bins_ttbar, dsigma_ttbar = ttbar_dsigma["delta_y_bins"]
plt.step(delta_y_bins_ttbar, dsigma_ttbar, where="mid", alpha=0.7,
         label=r"LHeC@1.2 TeV : $\gamma\gamma \to t\bar{t} (\times 10^{3})$", color="purple", linewidth=3)

# ✅ Plot τ⁺τ⁻ Background
delta_y_bins_tautau, dsigma_tautau = tautau_dsigma["delta_y_bins"]
plt.step(delta_y_bins_tautau, dsigma_tautau, where="mid", alpha=0.7,
         label=r"LHeC@1.2 TeV : $\gamma\gamma \to \tau^+\tau^-$", color="orange", linewidth=3)

# ✅ Plot μ⁺μ⁻ Background
delta_y_bins_mumu, dsigma_mumu = mumu_dsigma["delta_y_bins"]
plt.step(delta_y_bins_mumu, dsigma_mumu, where="mid", alpha=0.7,
         label=r"LHeC@1.2 TeV : $\gamma\gamma \to \mu^+\mu^-$", color="brown", linewidth=3)


# ✅ Inclusive tt̄ Background
delta_y_bins_ttbar_inc, dsigma_ttbar_inc = ttbar_inc_dsigma["delta_y_bins"]
plt.step(delta_y_bins_ttbar_inc, dsigma_ttbar_inc, where="mid", alpha=0.7,
         label=r"LHeC@1.2 TeV : Inclusive $t\bar{t}$", color="orchid", linewidth=3)

# ✅ Single Top Background
delta_y_bins_single_top, dsigma_single_top = single_top_dsigma["delta_y_bins"]
plt.step(delta_y_bins_single_top, dsigma_single_top, where="mid", alpha=0.7,
         label=r"LHeC@1.2 TeV : Single Top", color="teal", linewidth=3)

# ✅ W Production Background
delta_y_bins_w, dsigma_w = w_dsigma["delta_y_bins"]
plt.step(delta_y_bins_w, dsigma_w, where="mid", alpha=0.7,
         label=r"LHeC@1.2 TeV : W Production", color="darkgreen", linewidth=3)

# ✅ Z Production Background
delta_y_bins_z, dsigma_z = z_dsigma["delta_y_bins"]
plt.step(delta_y_bins_z, dsigma_z, where="mid", alpha=0.7,
         label=r"LHeC@1.2 TeV : Z Production", color="darkred", linewidth=3)


# Set labels and title
plt.xlabel(r"$\Delta y(\ell_1, \ell_2)$")
plt.ylabel(r"$\frac{d\sigma}{d\Delta y} \ \mathrm{[pb]}$")
plt.title(r"Delphes simulation : $e^- p \to e^- W^+ W^- p \to e^- \ell^+ \nu_{\ell} \ell^- \bar{\nu}_{\ell} p$ : LHeC@1.2 TeV", fontsize=20)
plt.yscale("log")

# Grid, legend, and layout adjustments
plt.legend()
plt.grid(True, linestyle="--", alpha=0.6)
plt.tight_layout()
plt.ylim(0.00001, 10.0)

# Save the plot
plt.savefig("/home/hamzeh-khanpour/Documents/GitHub/LHeC_Fast_Simulation/Pythia8_Delphes_fully_leptonic_allsignal_bkgs/differential_cross_section_deltaY_ll_allFMsignal_allbkgs.pdf", dpi=600)

# Show the plot
plt.show()




# ===================================================
# ✅ Missing Transverse Energy (MET) Differential Cross-Section
# ===================================================

plt.figure(figsize=(11, 12))  # Create a new figure for MET plot
plt.subplots_adjust(left=0.15, right=0.95, bottom=0.12, top=0.95)

# Loop through all signals and plot their differential cross-section for MET
for signal_name, dsigma_data in signal_dsigma.items():
    met_bins, dsigma = dsigma_data["met_bins"]
    plt.step(met_bins, dsigma, where="mid", alpha=0.7,
             label=f"LHeC@1.2 TeV : Signal ($W^+ W^-$) [{signal_name}]",
             color=signal_colors.get(signal_name, "black"), linewidth=3)

# ✅ Plot Background for MET
met_bins_background, dsigma_background = background_dsigma["met_bins"]
plt.step(met_bins_background, dsigma_background, where="mid", alpha=0.7,
         label=r"LHeC@1.2 TeV : SM ($\gamma\gamma \to W^+ W^-$)", color="blue", linewidth=3)


# ✅ Plot tt̄ Background
met_bins_ttbar, dsigma_ttbar = ttbar_dsigma["met_bins"]
plt.step(met_bins_ttbar, dsigma_ttbar, where="mid", alpha=0.7,
         label=r"LHeC@1.2 TeV : $\gamma\gamma \to t\bar{t} (\times 10^{3})$", color="purple", linewidth=3)


# ✅ τ⁺τ⁻ Background
met_bins_tautau, dsigma_tautau = tautau_dsigma["met_bins"]
plt.step(met_bins_tautau, dsigma_tautau, where="mid", alpha=0.7,
         label=r"LHeC@1.2 TeV : $\gamma\gamma \to \tau^+\tau^-$", color="orange", linewidth=3)

# ✅ μ⁺μ⁻ Background
met_bins_mumu, dsigma_mumu = mumu_dsigma["met_bins"]
plt.step(met_bins_mumu, dsigma_mumu, where="mid", alpha=0.7,
         label=r"LHeC@1.2 TeV : $\gamma\gamma \to \mu^+\mu^-$", color="brown", linewidth=3)

# ✅ Inclusive tt̄ Background
met_bins_ttbar_inc, dsigma_ttbar_inc = ttbar_inc_dsigma["met_bins"]
plt.step(met_bins_ttbar_inc, dsigma_ttbar_inc, where="mid", alpha=0.7,
         label=r"LHeC@1.2 TeV : Inclusive $t\bar{t}$", color="orchid", linewidth=3)

# ✅ Single Top Background
met_bins_single_top, dsigma_single_top = single_top_dsigma["met_bins"]
plt.step(met_bins_single_top, dsigma_single_top, where="mid", alpha=0.7,
         label=r"LHeC@1.2 TeV : Single Top", color="teal", linewidth=3)

# ✅ W Production Background
met_bins_w, dsigma_w = w_dsigma["met_bins"]
plt.step(met_bins_w, dsigma_w, where="mid", alpha=0.7,
         label=r"LHeC@1.2 TeV : W Production", color="darkgreen", linewidth=3)

# ✅ Z Production Background
met_bins_z, dsigma_z = z_dsigma["met_bins"]
plt.step(met_bins_z, dsigma_z, where="mid", alpha=0.7,
         label=r"LHeC@1.2 TeV : Z Production", color="darkred", linewidth=3)


# Set labels and title
plt.xlabel(r"$E_T^{miss} \ \mathrm{[GeV]}$")
plt.ylabel(r"$\frac{d\sigma}{dE_T^{miss}} \ \mathrm{[pb/GeV]}$")
plt.title(r"Delphes simulation : $e^- p \to e^- W^+ W^- p \to e^- \ell^+ \nu_{\ell} \ell^- \bar{\nu}_{\ell} p$ : LHeC@1.2 TeV", fontsize=20)
plt.yscale("log")

# Grid, legend, and layout adjustments
plt.legend()
plt.grid(True, linestyle="--", alpha=0.6)
plt.tight_layout()
plt.ylim(0.00001, 0.001)

# Save the plot
plt.savefig("/home/hamzeh-khanpour/Documents/GitHub/LHeC_Fast_Simulation/Pythia8_Delphes_fully_leptonic_allsignal_bkgs/differential_cross_section_MET_allFMsignal_allbkgs.pdf", dpi=600)

# Show the plot
plt.show()



# ===================================================
# ✅ Transverse Mass of W Boson \( M_T(W) \) Differential Cross-Section
# ===================================================

plt.figure(figsize=(11, 12))  # Create a new figure for M_T(W) plot
plt.subplots_adjust(left=0.15, right=0.95, bottom=0.12, top=0.95)

# Loop through all signals and plot their differential cross-section for M_T(W)
for signal_name, dsigma_data in signal_dsigma.items():
    mT_W_bins, dsigma = dsigma_data["mT_W_bins"]
    plt.step(mT_W_bins, dsigma, where="mid", alpha=0.7,
             label=f"LHeC@1.2 TeV : Signal ($W^+ W^-$) [{signal_name}]",
             color=signal_colors.get(signal_name, "black"), linewidth=3)

# ✅ Plot Background for M_T(W)
mT_W_bins_background, dsigma_background = background_dsigma["mT_W_bins"]
plt.step(mT_W_bins_background, dsigma_background, where="mid", alpha=0.7,
         label=r"LHeC@1.2 TeV : SM ($\gamma\gamma \to W^+ W^-$)", color="blue", linewidth=3)


# ✅ Plot tt̄ Background
mT_W_bins_ttbar, dsigma_ttbar = ttbar_dsigma["mT_W_bins"]
plt.step(mT_W_bins_ttbar, dsigma_ttbar, where="mid", alpha=0.7,
         label=r"LHeC@1.2 TeV : $\gamma\gamma \to t\bar{t} (\times 10^{3})$", color="purple", linewidth=3)



# ✅ τ⁺τ⁻ Background
mT_W_bins_tautau, dsigma_tautau = tautau_dsigma["mT_W_bins"]
plt.step(mT_W_bins_tautau, dsigma_tautau, where="mid", alpha=0.7,
         label=r"LHeC@1.2 TeV : $\gamma\gamma \to \tau^+\tau^-$", color="orange", linewidth=3)

# ✅ μ⁺μ⁻ Background
mT_W_bins_mumu, dsigma_mumu = mumu_dsigma["mT_W_bins"]
plt.step(mT_W_bins_mumu, dsigma_mumu, where="mid", alpha=0.7,
         label=r"LHeC@1.2 TeV : $\gamma\gamma \to \mu^+\mu^-$", color="brown", linewidth=3)

# ✅ Inclusive tt̄ Background
mT_W_bins_ttbar_inc, dsigma_ttbar_inc = ttbar_inc_dsigma["mT_W_bins"]
plt.step(mT_W_bins_ttbar_inc, dsigma_ttbar_inc, where="mid", alpha=0.7,
         label=r"LHeC@1.2 TeV : Inclusive $t\bar{t}$", color="orchid", linewidth=3)

# ✅ Single Top Background
mT_W_bins_single_top, dsigma_single_top = single_top_dsigma["mT_W_bins"]
plt.step(mT_W_bins_single_top, dsigma_single_top, where="mid", alpha=0.7,
         label=r"LHeC@1.2 TeV : Single Top", color="teal", linewidth=3)

# ✅ W Production Background
mT_W_bins_w, dsigma_w = w_dsigma["mT_W_bins"]
plt.step(mT_W_bins_w, dsigma_w, where="mid", alpha=0.7,
         label=r"LHeC@1.2 TeV : W Production", color="darkgreen", linewidth=3)

# ✅ Z Production Background
mT_W_bins_z, dsigma_z = z_dsigma["mT_W_bins"]
plt.step(mT_W_bins_z, dsigma_z, where="mid", alpha=0.7,
         label=r"LHeC@1.2 TeV : Z Production", color="darkred", linewidth=3)


# Set labels and title
plt.xlabel(r"$M_T(W) \ \mathrm{[GeV]}$")
plt.ylabel(r"$\frac{d\sigma}{dM_T(W)} \ \mathrm{[pb/GeV]}$")
plt.title(r"Delphes simulation : $e^- p \to e^- W^+ W^- p \to e^- \ell^+ \nu_{\ell} \ell^- \bar{\nu}_{\ell} p$ : LHeC@1.2 TeV", fontsize=20)
plt.yscale("log")

# Grid, legend, and layout adjustments
plt.legend()
plt.grid(True, linestyle="--", alpha=0.6)
plt.tight_layout()
plt.ylim(0.00001, 0.001)

# Save the plot
plt.savefig("/home/hamzeh-khanpour/Documents/GitHub/LHeC_Fast_Simulation/Pythia8_Delphes_fully_leptonic_allsignal_bkgs/differential_cross_section_mT_W_allFMsignal_allbkgs.pdf", dpi=600)

# Show the plot
plt.show()




# ===================================================
# ✅ Diboson Invariant Mass \( M_{WW} \) Differential Cross-Section
# ===================================================

plt.figure(figsize=(11, 12))  # Create a new figure for M_WW plot
plt.subplots_adjust(left=0.15, right=0.95, bottom=0.12, top=0.95)

# Loop through all signals and plot their differential cross-section for M_WW
for signal_name, dsigma_data in signal_dsigma.items():
    m_WW_bins, dsigma = dsigma_data["m_WW_bins"]
    plt.step(m_WW_bins, dsigma, where="mid", alpha=0.7,
             label=f"LHeC@1.2 TeV : Signal ($W^+ W^-$) [{signal_name}]",
             color=signal_colors.get(signal_name, "black"), linewidth=3)

# ✅ Plot Background for Diboson Invariant Mass
m_WW_bins_background, dsigma_background = background_dsigma["m_WW_bins"]
plt.step(m_WW_bins_background, dsigma_background, where="mid", alpha=0.7,
         label=r"LHeC@1.2 TeV : SM ($\gamma\gamma \to W^+ W^-$)", color="blue", linewidth=3)



# ✅ Plot tt̄ Background
m_WW_bins_ttbar, dsigma_ttbar = ttbar_dsigma["m_WW_bins"]
plt.step(m_WW_bins_ttbar, dsigma_ttbar, where="mid", alpha=0.7,
         label=r"LHeC@1.2 TeV : $\gamma\gamma \to t\bar{t} (\times 10^{3})$", color="purple", linewidth=3)


# ✅ Plot τ⁺τ⁻ Background
m_WW_bins_tautau, dsigma_tautau = tautau_dsigma["m_WW_bins"]
plt.step(m_WW_bins_tautau, dsigma_tautau, where="mid", alpha=0.7,
         label=r"LHeC@1.2 TeV : $\gamma\gamma \to \tau^+\tau^-$", color="orange", linewidth=3)

# ✅ Plot μ⁺μ⁻ Background
m_WW_bins_mumu, dsigma_mumu = mumu_dsigma["m_WW_bins"]
plt.step(m_WW_bins_mumu, dsigma_mumu, where="mid", alpha=0.7,
         label=r"LHeC@1.2 TeV : $\gamma\gamma \to \mu^+\mu^-$", color="brown", linewidth=3)


# ✅ Inclusive tt̄ Background
m_WW_bins_ttbar_inc, dsigma_ttbar_inc = ttbar_inc_dsigma["m_WW_bins"]
plt.step(m_WW_bins_ttbar_inc, dsigma_ttbar_inc, where="mid", alpha=0.7,
         label=r"LHeC@1.2 TeV : Inclusive $t\bar{t}$", color="orchid", linewidth=3)

# ✅ Single Top Background
m_WW_bins_single_top, dsigma_single_top = single_top_dsigma["m_WW_bins"]
plt.step(m_WW_bins_single_top, dsigma_single_top, where="mid", alpha=0.7,
         label=r"LHeC@1.2 TeV : Single Top", color="teal", linewidth=3)

# ✅ W Production Background
m_WW_bins_w, dsigma_w = w_dsigma["m_WW_bins"]
plt.step(m_WW_bins_w, dsigma_w, where="mid", alpha=0.7,
         label=r"LHeC@1.2 TeV : W Production", color="darkgreen", linewidth=3)

# ✅ Z Production Background
m_WW_bins_z, dsigma_z = z_dsigma["m_WW_bins"]
plt.step(m_WW_bins_z, dsigma_z, where="mid", alpha=0.7,
         label=r"LHeC@1.2 TeV : Z Production", color="darkred", linewidth=3)


# Set labels and title
plt.xlabel(r"$M_{WW} \ \mathrm{[GeV]}$")
plt.ylabel(r"$\frac{d\sigma}{dM_{WW}} \ \mathrm{[pb/GeV]}$")
plt.title(r"Delphes simulation : $e^- p \to e^- W^+ W^- p \to e^- \ell^+ \nu_{\ell} \ell^- \bar{\nu}_{\ell} p$ : LHeC@1.2 TeV", fontsize=20)
plt.yscale("log")

# Grid, legend, and layout adjustments
plt.legend()
plt.grid(True, linestyle="--", alpha=0.6)
plt.tight_layout()
plt.ylim(0.00001, 10.0)

# Save the plot
plt.savefig("/home/hamzeh-khanpour/Documents/GitHub/LHeC_Fast_Simulation/Pythia8_Delphes_fully_leptonic_allsignal_bkgs/differential_cross_section_m_WW_allFMsignal_allbkgs.pdf", dpi=600)

# Show the plot
plt.show()





# ===================================================
# ✅ Diboson Transverse Momentum \( p_T^{WW} \) Differential Cross-Section
# ===================================================

plt.figure(figsize=(11, 12))  # Create a new figure for pT_WW plot
plt.subplots_adjust(left=0.15, right=0.95, bottom=0.12, top=0.95)

# Loop through all signals and plot their differential cross-section for pT_WW
for signal_name, dsigma_data in signal_dsigma.items():
    pT_WW_bins, dsigma = dsigma_data["pT_WW_bins"]
    plt.step(pT_WW_bins, dsigma, where="mid", alpha=0.7,
             label=f"LHeC@1.2 TeV : Signal ($W^+ W^-$) [{signal_name}]",
             color=signal_colors.get(signal_name, "black"), linewidth=3)

# ✅ Plot Background for Diboson Transverse Momentum
pT_WW_bins_background, dsigma_background = background_dsigma["pT_WW_bins"]
plt.step(pT_WW_bins_background, dsigma_background, where="mid", alpha=0.7,
         label=r"LHeC@1.2 TeV : SM ($\gamma\gamma \to W^+ W^-$)", color="blue", linewidth=3)


# ✅ Plot tt̄ Background
pT_WW_bins_ttbar, dsigma_ttbar = ttbar_dsigma["pT_WW_bins"]
plt.step(pT_WW_bins_ttbar, dsigma_ttbar, where="mid", alpha=0.7,
         label=r"LHeC@1.2 TeV : $\gamma\gamma \to t\bar{t} (\times 10^{3})$", color="purple", linewidth=3)



# ✅ Plot τ⁺τ⁻ Background
pT_WW_bins_tautau, dsigma_tautau = tautau_dsigma["pT_WW_bins"]
plt.step(pT_WW_bins_tautau, dsigma_tautau, where="mid", alpha=0.7,
         label=r"LHeC@1.2 TeV : $\gamma\gamma \to \tau^+\tau^-$", color="orange", linewidth=3)

# ✅ Plot μ⁺μ⁻ Background
pT_WW_bins_mumu, dsigma_mumu = mumu_dsigma["pT_WW_bins"]
plt.step(pT_WW_bins_mumu, dsigma_mumu, where="mid", alpha=0.7,
         label=r"LHeC@1.2 TeV : $\gamma\gamma \to \mu^+\mu^-$", color="brown", linewidth=3)

# ✅ Inclusive tt̄ Background
pT_WW_bins_ttbar_inc, dsigma_ttbar_inc = ttbar_inc_dsigma["pT_WW_bins"]
plt.step(pT_WW_bins_ttbar_inc, dsigma_ttbar_inc, where="mid", alpha=0.7,
         label=r"LHeC@1.2 TeV : Inclusive $t\bar{t}$", color="orchid", linewidth=3)

# ✅ Single Top Background
pT_WW_bins_single_top, dsigma_single_top = single_top_dsigma["pT_WW_bins"]
plt.step(pT_WW_bins_single_top, dsigma_single_top, where="mid", alpha=0.7,
         label=r"LHeC@1.2 TeV : Single Top", color="teal", linewidth=3)

# ✅ W Production Background
pT_WW_bins_w, dsigma_w = w_dsigma["pT_WW_bins"]
plt.step(pT_WW_bins_w, dsigma_w, where="mid", alpha=0.7,
         label=r"LHeC@1.2 TeV : W Production", color="darkgreen", linewidth=3)

# ✅ Z Production Background
pT_WW_bins_z, dsigma_z = z_dsigma["pT_WW_bins"]
plt.step(pT_WW_bins_z, dsigma_z, where="mid", alpha=0.7,
         label=r"LHeC@1.2 TeV : Z Production", color="darkred", linewidth=3)


# Set labels and title
plt.xlabel(r"$p_T^{WW} \ \mathrm{[GeV]}$")
plt.ylabel(r"$\frac{d\sigma}{dp_T^{WW}} \ \mathrm{[pb/GeV]}$")
plt.title(r"Delphes simulation : $e^- p \to e^- W^+ W^- p \to e^- \ell^+ \nu_{\ell} \ell^- \bar{\nu}_{\ell} p$ : LHeC@1.2 TeV", fontsize=20)
plt.yscale("log")

# Grid, legend, and layout adjustments
plt.legend()
plt.grid(True, linestyle="--", alpha=0.6)
plt.tight_layout()
plt.ylim(0.00001, 0.001)

# Save the plot
plt.savefig("/home/hamzeh-khanpour/Documents/GitHub/LHeC_Fast_Simulation/Pythia8_Delphes_fully_leptonic_allsignal_bkgs/differential_cross_section_pT_WW_allFMsignal_allbkgs.pdf", dpi=600)

# Show the plot
plt.show()




# ===================================================
# ✅ Leading Lepton Pseudorapidity \( \eta_{\mathrm{leading}~\ell} \) Differential Cross-Section
# ===================================================

plt.figure(figsize=(11, 12))  # Create a new figure for leading lepton eta plot
plt.subplots_adjust(left=0.15, right=0.95, bottom=0.12, top=0.95)

# Loop through all signals and plot their differential cross-section for leading lepton eta
for signal_name, dsigma_data in signal_dsigma.items():
    eta_leading_bins, dsigma = dsigma_data["eta_bins_leading_lepton"]
    plt.step(eta_leading_bins, dsigma, where="mid", alpha=0.7,
             label=f"LHeC@1.2 TeV : Signal ($W^+ W^-$) [{signal_name}]",
             color=signal_colors.get(signal_name, "black"), linewidth=3)

# ✅ Plot Background for Leading Lepton Eta
eta_leading_bins_background, dsigma_background = background_dsigma["eta_bins_leading_lepton"]
plt.step(eta_leading_bins_background, dsigma_background, where="mid", alpha=0.7,
         label=r"LHeC@1.2 TeV : SM ($\gamma\gamma \to W^+ W^-$)", color="blue", linewidth=3)


# ✅ Plot tt̄ Background
eta_leading_bins_ttbar, dsigma_ttbar = ttbar_dsigma["eta_bins_leading_lepton"]
plt.step(eta_leading_bins_ttbar, dsigma_ttbar, where="mid", alpha=0.7,
         label=r"LHeC@1.2 TeV : $\gamma\gamma \to t\bar{t} (\times 10^{3})$", color="purple", linewidth=3)



# ✅ τ⁺τ⁻ Background
eta_leading_bins_tautau, dsigma_tautau = tautau_dsigma["eta_bins_leading_lepton"]
plt.step(eta_leading_bins_tautau, dsigma_tautau, where="mid", alpha=0.7,
         label=r"LHeC@1.2 TeV : $\gamma\gamma \to \tau^+\tau^-$", color="orange", linewidth=3)

# ✅ μ⁺μ⁻ Background
eta_leading_bins_mumu, dsigma_mumu = mumu_dsigma["eta_bins_leading_lepton"]
plt.step(eta_leading_bins_mumu, dsigma_mumu, where="mid", alpha=0.7,
         label=r"LHeC@1.2 TeV : $\gamma\gamma \to \mu^+\mu^-$", color="brown", linewidth=3)

# ✅ Inclusive tt̄ Background
eta_leading_bins_ttbar_inc, dsigma_ttbar_inc = ttbar_inc_dsigma["eta_bins_leading_lepton"]
plt.step(eta_leading_bins_ttbar_inc, dsigma_ttbar_inc, where="mid", alpha=0.7,
         label=r"LHeC@1.2 TeV : Inclusive $t\bar{t}$", color="orchid", linewidth=3)

# ✅ Single Top Background
eta_leading_bins_single_top, dsigma_single_top = single_top_dsigma["eta_bins_leading_lepton"]
plt.step(eta_leading_bins_single_top, dsigma_single_top, where="mid", alpha=0.7,
         label=r"LHeC@1.2 TeV : Single Top", color="teal", linewidth=3)

# ✅ W Production Background
eta_leading_bins_w, dsigma_w = w_dsigma["eta_bins_leading_lepton"]
plt.step(eta_leading_bins_w, dsigma_w, where="mid", alpha=0.7,
         label=r"LHeC@1.2 TeV : W Production", color="darkgreen", linewidth=3)

# ✅ Z Production Background
eta_leading_bins_z, dsigma_z = z_dsigma["eta_bins_leading_lepton"]
plt.step(eta_leading_bins_z, dsigma_z, where="mid", alpha=0.7,
         label=r"LHeC@1.2 TeV : Z Production", color="darkred", linewidth=3)

# Set labels and title
plt.xlabel(r"$\eta_{\mathrm{leading}~\ell}$")
plt.ylabel(r"$\frac{d\sigma}{d\eta_{\mathrm{leading}~\ell}} \ \mathrm{[pb]}$")
plt.title(r"Delphes simulation : $e^- p \to e^- W^+ W^- p \to e^- \ell^+ \nu_{\ell} \ell^- \bar{\nu}_{\ell} p$ : LHeC@1.2 TeV", fontsize=20)

# Grid, legend, and layout adjustments
plt.legend()
plt.grid(True, linestyle="--", alpha=0.6)
plt.tight_layout()
plt.ylim(0.0005, 0.003)

# Save the plot
plt.savefig("/home/hamzeh-khanpour/Documents/GitHub/LHeC_Fast_Simulation/Pythia8_Delphes_fully_leptonic_allsignal_bkgs/differential_cross_section_eta_leading_lepton_allFMsignal_allbkgs.pdf", dpi=600)

# Show the plot
plt.show()





# ===================================================
# ✅ Subleading Lepton Pseudorapidity \( \eta_^{\mathrm{subleading}~\ell} \) Differential Cross-Section
# ===================================================

plt.figure(figsize=(11, 12))  # Create a new figure for subleading lepton eta plot
plt.subplots_adjust(left=0.15, right=0.95, bottom=0.12, top=0.95)

# Loop through all signals and plot their differential cross-section for subleading lepton eta
for signal_name, dsigma_data in signal_dsigma.items():
    eta_subleading_bins, dsigma = dsigma_data["eta_bins_subleading_lepton"]
    plt.step(eta_subleading_bins, dsigma, where="mid", alpha=0.7,
             label=f"LHeC@1.2 TeV : Signal ($W^+ W^-$) [{signal_name}]",
             color=signal_colors.get(signal_name, "black"), linewidth=3)

# ✅ Plot Background for Subleading Lepton Eta
eta_subleading_bins_background, dsigma_background = background_dsigma["eta_bins_subleading_lepton"]
plt.step(eta_subleading_bins_background, dsigma_background, where="mid", alpha=0.7,
         label=r"LHeC@1.2 TeV : SM ($\gamma\gamma \to W^+ W^-$)", color="blue", linewidth=3)


# ✅ Plot tt̄ Background
eta_subleading_bins_ttbar, dsigma_ttbar = ttbar_dsigma["eta_bins_subleading_lepton"]
plt.step(eta_subleading_bins_ttbar, dsigma_ttbar, where="mid", alpha=0.7,
         label=r"LHeC@1.2 TeV : $\gamma\gamma \to t\bar{t} (\times 10^{3})$", color="purple", linewidth=3)


# ✅ Plot τ⁺τ⁻ Background
eta_subleading_bins_tautau, dsigma_tautau = tautau_dsigma["eta_bins_subleading_lepton"]
plt.step(eta_subleading_bins_tautau, dsigma_tautau, where="mid", alpha=0.7,
         label=r"LHeC@1.2 TeV : $\gamma\gamma \to \tau^+\tau^-$", color="orange", linewidth=3)

# ✅ Plot μ⁺μ⁻ Background
eta_subleading_bins_mumu, dsigma_mumu = mumu_dsigma["eta_bins_subleading_lepton"]
plt.step(eta_subleading_bins_mumu, dsigma_mumu, where="mid", alpha=0.7,
         label=r"LHeC@1.2 TeV : $\gamma\gamma \to \mu^+\mu^-$", color="brown", linewidth=3)


# ✅ Inclusive tt̄ Background
eta_subleading_bins_ttbar_inc, dsigma_ttbar_inc = ttbar_inc_dsigma["eta_bins_subleading_lepton"]
plt.step(eta_subleading_bins_ttbar_inc, dsigma_ttbar_inc, where="mid", alpha=0.7,
         label=r"LHeC@1.2 TeV : Inclusive $t\bar{t}$", color="orchid", linewidth=3)

# ✅ Single Top Background
eta_subleading_bins_single_top, dsigma_single_top = single_top_dsigma["eta_bins_subleading_lepton"]
plt.step(eta_subleading_bins_single_top, dsigma_single_top, where="mid", alpha=0.7,
         label=r"LHeC@1.2 TeV : Single Top", color="teal", linewidth=3)

# ✅ W Production Background
eta_subleading_bins_w, dsigma_w = w_dsigma["eta_bins_subleading_lepton"]
plt.step(eta_subleading_bins_w, dsigma_w, where="mid", alpha=0.7,
         label=r"LHeC@1.2 TeV : W Production", color="darkgreen", linewidth=3)

# ✅ Z Production Background
eta_subleading_bins_z, dsigma_z = z_dsigma["eta_bins_subleading_lepton"]
plt.step(eta_subleading_bins_z, dsigma_z, where="mid", alpha=0.7,
         label=r"LHeC@1.2 TeV : Z Production", color="darkred", linewidth=3)


# Set labels and title
plt.xlabel(r"$\eta_{\mathrm{subleading}~\ell}$")
plt.ylabel(r"$\frac{d\sigma}{d\eta_{\mathrm{subleading}~\ell}} \ \mathrm{[pb]}$")
plt.title(r"Delphes simulation : $e^- p \to e^- W^+ W^- p \to e^- \ell^+ \nu_{\ell} \ell^- \bar{\nu}_{\ell} p$ : LHeC@1.2 TeV", fontsize=20)

# Grid, legend, and layout adjustments
plt.legend()
plt.grid(True, linestyle="--", alpha=0.6)
plt.tight_layout()
plt.ylim(0.0005, 0.003)


# Save the plot
plt.savefig("/home/hamzeh-khanpour/Documents/GitHub/LHeC_Fast_Simulation/Pythia8_Delphes_fully_leptonic_allsignal_bkgs/differential_cross_section_eta_subleading_lepton_allFMsignal_allbkgs.pdf", dpi=600)

# Show the plot
plt.show()






#=========================================================================
#=========================================================================

# Open the output ROOT file (use "RECREATE" to overwrite, or "UPDATE" to append)
#output_file = ROOT.TFile("/home/hamzeh-khanpour/Documents/GitHub/LHeC_Fast_Simulation/Pythia8_Delphes_fully_leptonic_allsignal_bkgs/output_histograms_fully_leptonic_allbkgs.root", "RECREATE")

try:
    # ✅ Fix 1: Rename signal names to be ROOT-compatible (remove LaTeX)
    signal_dirs = {}
    for signal_name in signal_files.keys():
        clean_signal_name = signal_name.replace("$", "").replace("{", "").replace("}", "").replace("\\", "").replace(" ", "_")
        signal_dirs[signal_name] = output_file.mkdir(f"signal_{clean_signal_name}")

    # ✅ Fix 2: Create Background Directory BEFORE Writing Any Data
    background_dir = output_file.mkdir("SM_background")

    # ✅ Fix 3: Ensure `background_histograms` is properly defined for Fully Leptonic Analysis
    background_histograms = {
        "hist_leading_lepton": hist_leading_lepton_background,
        "hist_subleading_lepton": hist_subleading_lepton_background,
        "hist_m_ll": hist_m_ll_background,
        "hist_pt_ll": hist_pt_ll_background,
        "hist_rapidity_ll": hist_rapidity_ll_background,
        "hist_lepton_separation": hist_lepton_separation_background,
        "hist_azimuthal_angle_separation": hist_azimuthal_angle_separation_background,
        "hist_rapidity_difference": hist_rapidity_difference_background,
        "hist_met": hist_met_background,
        "hist_mT_W": hist_mT_W_background,
        "hist_m_WW": hist_m_WW_background,
        "hist_pT_WW": hist_pT_WW_background,
        "hist_leading_lepton_eta": hist_leading_lepton_eta_background,
        "hist_subleading_lepton_eta": hist_subleading_lepton_eta_background
    }

    # ✅ Fix 4: Save signal histograms (DO NOT redefine `signal_histograms`)
    for signal_name, histograms in signal_histograms.items():
        if signal_name in signal_dirs:  # Ensure directory exists
            signal_dirs[signal_name].cd()
            for hist_name, hist in histograms.items():
                if hist:  # Ensure histogram exists
                    hist.Write()
                else:
                    print(f"⚠️ Warning: Histogram {hist_name} for {signal_name} is empty!")

    # ✅ Fix 5: Save background histograms
    background_dir.cd()
    for hist_name, hist in background_histograms.items():
        if hist:  # Ensure histogram exists
            hist.Write()
        else:
            print(f"⚠️ Warning: Background histogram {hist_name} is empty!")

    print("✅ Histograms successfully saved to output_histograms_fully_leptonic.root with separate branches for each signal and 'SM_background'.")

except Exception as e:
    print(f"❌ Error while saving histograms: {e}")

finally:
    # ✅ Fix 6: Always close the ROOT file properly
    output_file.Close()
    print("📁 ROOT file closed successfully.")

#=========================================================================
#=========================================================================
