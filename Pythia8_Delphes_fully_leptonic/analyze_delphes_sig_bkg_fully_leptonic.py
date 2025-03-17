
## ================================================================================
##        Hamzeh Khanpour  --- January 2025
## ================================================================================




#!/usr/bin/env python

import ROOT
from ROOT import TLorentzVector
import matplotlib.pyplot as plt
import numpy as np  # Import numpy for calculations


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
signal_file_path = "aa_ww_fully_leptonic_NP_1_FM2_100.root"
background_file_path = "aa_ww_fully_leptonic_SM.root"



# Load Delphes library
ROOT.gSystem.Load("libDelphes")
ROOT.gInterpreter.Declare('#include "classes/DelphesClasses.h"')
ROOT.gInterpreter.Declare('#include "external/ExRootAnalysis/ExRootTreeReader.h"')


#==============================================================================
#==============================================================================

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
    numberOfEntries = treeReader.GetEntries()

    # Counters for efficiency calculation
    total_events = numberOfEntries
    selected_events = 0

    # Get branches for electrons, muons, and MET
    branchElectron = treeReader.UseBranch("Electron")
    branchMuon = treeReader.UseBranch("Muon")
    branchMissingET = treeReader.UseBranch("MissingET")
    branchJet = treeReader.UseBranch("Jet")

    # Process each event
    for entry in range(numberOfEntries):
        treeReader.ReadEntry(entry)

        # Collect leptons (electrons + muons)
        leptons = []
        for i in range(branchElectron.GetEntries()):
            electron = branchElectron.At(i)
            lepton_vec = ROOT.TLorentzVector()
            lepton_vec.SetPtEtaPhiM(electron.PT, electron.Eta, electron.Phi, 0.0)
            if electron.PT > 10:
                leptons.append(lepton_vec)
        for i in range(branchMuon.GetEntries()):
            muon = branchMuon.At(i)
            lepton_vec = ROOT.TLorentzVector()
            lepton_vec.SetPtEtaPhiM(muon.PT, muon.Eta, muon.Phi, 0.0)
            if muon.PT > 10:
                leptons.append(lepton_vec)

        # Collect jets (for vetoing)
        jets = []
        for i in range(branchJet.GetEntries()):
            jet = branchJet.At(i)
            jet_vec = ROOT.TLorentzVector()
            jet_vec.SetPtEtaPhiM(jet.PT, jet.Eta, jet.Phi, jet.Mass)
            if jet.PT > 10:
                jets.append(jet_vec)

        # Apply selection criteria: exactly two leptons and no jets
        if len(leptons) != 2 or len(jets) > 0:
            continue

        # Count selected events
        selected_events += 1

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

        hist_lepton_separation.Fill(delta_r)  # ΔR(ℓ1, ℓ2)
        hist_azimuthal_angle_separation.Fill(delta_phi)  # Δφ(ℓ1, ℓ2)
        hist_rapidity_difference.Fill(delta_y)  # Δy(ℓ1, ℓ2)

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


    # Calculate selection efficiency
    efficiency = selected_events / total_events if total_events > 0 else 0


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
        efficiency,
    )




#==============================================================================
#==============================================================================



# Create histograms for fully leptonic kinematic distributions (signal and background)

# ✅ Leading Lepton pT Distribution
hist_leading_lepton_signal = ROOT.TH1F(
    "hist_leading_lepton_signal", "Leading Lepton pT Distribution; p_{T}^{leading} [GeV]; Entries", 50, 0, 400
)
hist_leading_lepton_background = ROOT.TH1F(
    "hist_leading_lepton_background", "Leading Lepton pT Distribution; p_{T}^{leading} [GeV]; Entries", 50, 0, 400
)

# ✅ Subleading Lepton pT Distribution
hist_subleading_lepton_signal = ROOT.TH1F(
    "hist_subleading_lepton_signal", "Subleading Lepton pT Distribution; p_{T}^{subleading} [GeV]; Entries", 50, 0, 300
)
hist_subleading_lepton_background = ROOT.TH1F(
    "hist_subleading_lepton_background", "Subleading Lepton pT Distribution; p_{T}^{subleading} [GeV]; Entries", 50, 0, 300
)

# ✅ Dilepton Invariant Mass Distribution (M_ll)
hist_m_ll_signal = ROOT.TH1F(
    "hist_m_ll_signal", "Dilepton Invariant Mass Distribution; M_{\\ell^+\\ell^-} [GeV]; Entries", 50, 0, 800
)
hist_m_ll_background = ROOT.TH1F(
    "hist_m_ll_background", "Dilepton Invariant Mass Distribution; M_{\\ell^+\\ell^-} [GeV]; Entries", 50, 0, 800
)

# ✅ Dilepton Transverse Momentum Distribution (pT_ll)
hist_pt_ll_signal = ROOT.TH1F(
    "hist_pt_ll_signal", "Dilepton Transverse Momentum Distribution; p_{T}^{\\ell^+\\ell^-} [GeV]; Entries", 50, 0, 300
)
hist_pt_ll_background = ROOT.TH1F(
    "hist_pt_ll_background", "Dilepton Transverse Momentum Distribution; p_{T}^{\\ell^+\\ell^-} [GeV]; Entries", 50, 0, 300
)

# ✅ Rapidity of Dilepton System (Y_ll)
hist_rapidity_ll_signal = ROOT.TH1F(
    "hist_rapidity_ll_signal", "Rapidity of Dilepton System; Y_{\\ell^+\\ell^-}; Entries", 50, -5, 5
)
hist_rapidity_ll_background = ROOT.TH1F(
    "hist_rapidity_ll_background", "Rapidity of Dilepton System; Y_{\\ell^+\\ell^-}; Entries", 50, -5, 5
)

# ✅ Leading Lepton Pseudorapidity (η)
hist_leading_lepton_eta_signal = ROOT.TH1F(
    "hist_leading_lepton_eta_signal", "Leading Lepton Pseudorapidity; \\eta^{leading}_{\\ell}; Entries", 50, -5, 5
)
hist_leading_lepton_eta_background = ROOT.TH1F(
    "hist_leading_lepton_eta_background", "Leading Lepton Pseudorapidity; \\eta^{leading}_{\\ell}; Entries", 50, -5, 5
)

# ✅ Subleading Lepton Pseudorapidity (η)
hist_subleading_lepton_eta_signal = ROOT.TH1F(
    "hist_subleading_lepton_eta_signal", "Subleading Lepton Pseudorapidity; \\eta^{subleading}_{\\ell}; Entries", 50, -5, 5
)
hist_subleading_lepton_eta_background = ROOT.TH1F(
    "hist_subleading_lepton_eta_background", "Subleading Lepton Pseudorapidity; \\eta^{subleading}_{\\ell}; Entries", 50, -5, 5
)





# ✅ Lepton Separation (ΔR)
hist_lepton_separation_signal = ROOT.TH1F(
    "hist_lepton_separation_signal", "Lepton Separation; \\Delta R(\\ell_1, \\ell_2); Entries", 50, 0, 5
)
hist_lepton_separation_background = ROOT.TH1F(
    "hist_lepton_separation_background", "Lepton Separation; \\Delta R(\\ell_1, \\ell_2); Entries", 50, 0, 5
)

# ✅ Azimuthal Angle Separation (Δφ)
hist_azimuthal_angle_separation_signal = ROOT.TH1F(
    "hist_azimuthal_angle_separation_signal", "Azimuthal Angle Separation; \\Delta \\phi(\\ell_1, \\ell_2); Entries", 50, 0, 3.2
)
hist_azimuthal_angle_separation_background = ROOT.TH1F(
    "hist_azimuthal_angle_separation_background", "Azimuthal Angle Separation; \\Delta \\phi(\\ell_1, \\ell_2); Entries", 50, 0, 3.2
)

# ✅ Rapidity Difference (Δy)
hist_rapidity_difference_signal = ROOT.TH1F(
    "hist_rapidity_difference_signal", "Rapidity Difference; \\Delta y(\\ell_1, \\ell_2); Entries", 50, 0, 5
)
hist_rapidity_difference_background = ROOT.TH1F(
    "hist_rapidity_difference_background", "Rapidity Difference; \\Delta y(\\ell_1, \\ell_2); Entries", 50, 0, 5
)


# ✅ Missing Transverse Energy (MET)
hist_met_signal = ROOT.TH1F(
    "hist_met_signal", "Missing Transverse Energy; MET [GeV]; Entries", 50, 0, 400
)
hist_met_background = ROOT.TH1F(
    "hist_met_background", "Missing Transverse Energy; MET [GeV]; Entries", 50, 0, 400
)




# ✅ Transverse Mass of W boson (M_T(W))
hist_mT_W_signal = ROOT.TH1F(
    "hist_mT_W_signal", "Transverse Mass of W Boson; M_T(W) [GeV]; Entries", 50, 0, 500
)
hist_mT_W_background = ROOT.TH1F(
    "hist_mT_W_background", "Transverse Mass of W Boson; M_T(W) [GeV]; Entries", 50, 0, 00
)

# ✅ Reconstructed Diboson System Mass (M_WW)
hist_m_WW_signal = ROOT.TH1F(
    "hist_m_WW_signal", "Reconstructed Diboson Mass; M_{WW} [GeV]; Entries", 50, 0, 1000
)
hist_m_WW_background = ROOT.TH1F(
    "hist_m_WW_background", "Reconstructed Diboson Mass; M_{WW} [GeV]; Entries", 50, 0, 1000
)

# ✅ Diboson Transverse Momentum (pT_WW)
hist_pT_WW_signal = ROOT.TH1F(
    "hist_pT_WW_signal", "Diboson Transverse Momentum; p_T^{WW} [GeV]; Entries", 50, 0, 400
)
hist_pT_WW_background = ROOT.TH1F(
    "hist_pT_WW_background", "Diboson Transverse Momentum; p_T^{WW} [GeV]; Entries", 50, 0, 400
)







#==============================================================================
#==============================================================================




# Process signal and background files and calculate efficiencies
# Call the process_file function and update the unpacking to include all returned values

(hist_leading_lepton_signal,
 hist_subleading_lepton_signal,
 hist_m_ll_signal,
 hist_pt_ll_signal,
 hist_rapidity_ll_signal,
 hist_lepton_separation_signal,  # ✅ Lepton Separation ΔR (Signal)
 hist_azimuthal_angle_separation_signal,  # ✅ Azimuthal Angle Separation Δφ (Signal)
 hist_rapidity_difference_signal,  # ✅ Rapidity Difference Δy (Signal)
 hist_met_signal,  # ✅ Missing Transverse Energy (MET) (Signal)
 hist_mT_W_signal,  # ✅ Transverse Mass of W boson (Signal)
 hist_m_WW_signal,  # ✅ Reconstructed Diboson System Mass (Signal)
 hist_pT_WW_signal,  # ✅ Diboson Transverse Momentum (Signal)
 hist_leading_lepton_eta_signal,   # ✅ Leading Lepton η histogram (Signal)
 hist_subleading_lepton_eta_signal,  # ✅ Subleading Lepton η histogram (Signal)
 signal_efficiency) = process_file(
    signal_file_path,
    hist_leading_lepton_signal,
    hist_subleading_lepton_signal,
    hist_m_ll_signal,
    hist_pt_ll_signal,
    hist_rapidity_ll_signal,
    hist_lepton_separation_signal,  # ✅ Pass to function
    hist_azimuthal_angle_separation_signal,  # ✅ Pass to function
    hist_rapidity_difference_signal,  # ✅ Pass to function
    hist_met_signal,  # ✅ Pass to function
    hist_mT_W_signal,  # ✅ Pass to function
    hist_m_WW_signal,  # ✅ Pass to function
    hist_pT_WW_signal,  # ✅ Pass to function
    hist_leading_lepton_eta_signal,
    hist_subleading_lepton_eta_signal
)


(hist_leading_lepton_background,
 hist_subleading_lepton_background,
 hist_m_ll_background,
 hist_pt_ll_background,
 hist_rapidity_ll_background,
 hist_lepton_separation_background,  # ✅ Lepton Separation ΔR (Background)
 hist_azimuthal_angle_separation_background,  # ✅ Azimuthal Angle Separation Δφ (Background)
 hist_rapidity_difference_background,  # ✅ Rapidity Difference Δy (Background)
 hist_met_background,  # ✅ Missing Transverse Energy (MET) (Background)
 hist_mT_W_background,  # ✅ Transverse Mass of W boson (Background)
 hist_m_WW_background,  # ✅ Reconstructed Diboson System Mass (Background)
 hist_pT_WW_background,  # ✅ Diboson Transverse Momentum (Background)
 hist_leading_lepton_eta_background,   # ✅ Leading Lepton η histogram (Background)
 hist_subleading_lepton_eta_background,  # ✅ Subleading Lepton η histogram (Background)
 background_efficiency) = process_file(
    background_file_path,
    hist_leading_lepton_background,
    hist_subleading_lepton_background,
    hist_m_ll_background,
    hist_pt_ll_background,
    hist_rapidity_ll_background,
    hist_lepton_separation_background,  # ✅ Pass to function
    hist_azimuthal_angle_separation_background,  # ✅ Pass to function
    hist_rapidity_difference_background,  # ✅ Pass to function
    hist_met_background,  # ✅ Pass to function
    hist_mT_W_background,  # ✅ Pass to function
    hist_m_WW_background,  # ✅ Pass to function
    hist_pT_WW_background,  # ✅ Pass to function
    hist_leading_lepton_eta_background,
    hist_subleading_lepton_eta_background
)








# Print selection efficiencies
print(f"Signal Selection Efficiency: {signal_efficiency:.2%}")
print(f"Background Selection Efficiency: {background_efficiency:.2%}")




#==============================================================================
#==============================================================================



# Convert ROOT histograms to matplotlib for leading lepton pT
x_vals_leading_lepton_signal = [hist_leading_lepton_signal.GetBinCenter(i) for i in range(1, hist_leading_lepton_signal.GetNbinsX() + 1)]
y_vals_leading_lepton_signal = [hist_leading_lepton_signal.GetBinContent(i) for i in range(1, hist_leading_lepton_signal.GetNbinsX() + 1)]

x_vals_leading_lepton_background = [hist_leading_lepton_background.GetBinCenter(i) for i in range(1, hist_leading_lepton_background.GetNbinsX() + 1)]
y_vals_leading_lepton_background = [hist_leading_lepton_background.GetBinContent(i) for i in range(1, hist_leading_lepton_background.GetNbinsX() + 1)]

# Convert ROOT histograms to matplotlib for subleading lepton pT
x_vals_subleading_lepton_signal = [hist_subleading_lepton_signal.GetBinCenter(i) for i in range(1, hist_subleading_lepton_signal.GetNbinsX() + 1)]
y_vals_subleading_lepton_signal = [hist_subleading_lepton_signal.GetBinContent(i) for i in range(1, hist_subleading_lepton_signal.GetNbinsX() + 1)]

x_vals_subleading_lepton_background = [hist_subleading_lepton_background.GetBinCenter(i) for i in range(1, hist_subleading_lepton_background.GetNbinsX() + 1)]
y_vals_subleading_lepton_background = [hist_subleading_lepton_background.GetBinContent(i) for i in range(1, hist_subleading_lepton_background.GetNbinsX() + 1)]

# Convert ROOT histograms to matplotlib for dilepton invariant mass
x_vals_m_ll_signal = [hist_m_ll_signal.GetBinCenter(i) for i in range(1, hist_m_ll_signal.GetNbinsX() + 1)]
y_vals_m_ll_signal = [hist_m_ll_signal.GetBinContent(i) for i in range(1, hist_m_ll_signal.GetNbinsX() + 1)]

x_vals_m_ll_background = [hist_m_ll_background.GetBinCenter(i) for i in range(1, hist_m_ll_background.GetNbinsX() + 1)]
y_vals_m_ll_background = [hist_m_ll_background.GetBinContent(i) for i in range(1, hist_m_ll_background.GetNbinsX() + 1)]

# Convert ROOT histograms to matplotlib for dilepton transverse momentum
x_vals_pt_ll_signal = [hist_pt_ll_signal.GetBinCenter(i) for i in range(1, hist_pt_ll_signal.GetNbinsX() + 1)]
y_vals_pt_ll_signal = [hist_pt_ll_signal.GetBinContent(i) for i in range(1, hist_pt_ll_signal.GetNbinsX() + 1)]

x_vals_pt_ll_background = [hist_pt_ll_background.GetBinCenter(i) for i in range(1, hist_pt_ll_background.GetNbinsX() + 1)]
y_vals_pt_ll_background = [hist_pt_ll_background.GetBinContent(i) for i in range(1, hist_pt_ll_background.GetNbinsX() + 1)]

# Convert ROOT histograms to matplotlib for dilepton rapidity
x_vals_rapidity_ll_signal = [hist_rapidity_ll_signal.GetBinCenter(i) for i in range(1, hist_rapidity_ll_signal.GetNbinsX() + 1)]
y_vals_rapidity_ll_signal = [hist_rapidity_ll_signal.GetBinContent(i) for i in range(1, hist_rapidity_ll_signal.GetNbinsX() + 1)]

x_vals_rapidity_ll_background = [hist_rapidity_ll_background.GetBinCenter(i) for i in range(1, hist_rapidity_ll_background.GetNbinsX() + 1)]
y_vals_rapidity_ll_background = [hist_rapidity_ll_background.GetBinContent(i) for i in range(1, hist_rapidity_ll_background.GetNbinsX() + 1)]

# ✅ NEW: Convert ROOT histograms for leading lepton pseudorapidity (η)
x_vals_leading_lepton_eta_signal = [hist_leading_lepton_eta_signal.GetBinCenter(i) for i in range(1, hist_leading_lepton_eta_signal.GetNbinsX() + 1)]
y_vals_leading_lepton_eta_signal = [hist_leading_lepton_eta_signal.GetBinContent(i) for i in range(1, hist_leading_lepton_eta_signal.GetNbinsX() + 1)]

x_vals_leading_lepton_eta_background = [hist_leading_lepton_eta_background.GetBinCenter(i) for i in range(1, hist_leading_lepton_eta_background.GetNbinsX() + 1)]
y_vals_leading_lepton_eta_background = [hist_leading_lepton_eta_background.GetBinContent(i) for i in range(1, hist_leading_lepton_eta_background.GetNbinsX() + 1)]

# ✅ NEW: Convert ROOT histograms for subleading lepton pseudorapidity (η)
x_vals_subleading_lepton_eta_signal = [hist_subleading_lepton_eta_signal.GetBinCenter(i) for i in range(1, hist_subleading_lepton_eta_signal.GetNbinsX() + 1)]
y_vals_subleading_lepton_eta_signal = [hist_subleading_lepton_eta_signal.GetBinContent(i) for i in range(1, hist_subleading_lepton_eta_signal.GetNbinsX() + 1)]

x_vals_subleading_lepton_eta_background = [hist_subleading_lepton_eta_background.GetBinCenter(i) for i in range(1, hist_subleading_lepton_eta_background.GetNbinsX() + 1)]
y_vals_subleading_lepton_eta_background = [hist_subleading_lepton_eta_background.GetBinContent(i) for i in range(1, hist_subleading_lepton_eta_background.GetNbinsX() + 1)]



# ✅ Convert ROOT histograms for lepton separation (ΔR)
x_vals_lepton_separation_signal = [hist_lepton_separation_signal.GetBinCenter(i) for i in range(1, hist_lepton_separation_signal.GetNbinsX() + 1)]
y_vals_lepton_separation_signal = [hist_lepton_separation_signal.GetBinContent(i) for i in range(1, hist_lepton_separation_signal.GetNbinsX() + 1)]

x_vals_lepton_separation_background = [hist_lepton_separation_background.GetBinCenter(i) for i in range(1, hist_lepton_separation_background.GetNbinsX() + 1)]
y_vals_lepton_separation_background = [hist_lepton_separation_background.GetBinContent(i) for i in range(1, hist_lepton_separation_background.GetNbinsX() + 1)]

# ✅ Convert ROOT histograms for azimuthal angle separation (Δφ)
x_vals_azimuthal_angle_separation_signal = [hist_azimuthal_angle_separation_signal.GetBinCenter(i) for i in range(1, hist_azimuthal_angle_separation_signal.GetNbinsX() + 1)]
y_vals_azimuthal_angle_separation_signal = [hist_azimuthal_angle_separation_signal.GetBinContent(i) for i in range(1, hist_azimuthal_angle_separation_signal.GetNbinsX() + 1)]

x_vals_azimuthal_angle_separation_background = [hist_azimuthal_angle_separation_background.GetBinCenter(i) for i in range(1, hist_azimuthal_angle_separation_background.GetNbinsX() + 1)]
y_vals_azimuthal_angle_separation_background = [hist_azimuthal_angle_separation_background.GetBinContent(i) for i in range(1, hist_azimuthal_angle_separation_background.GetNbinsX() + 1)]

# ✅ Convert ROOT histograms for rapidity difference (Δy)
x_vals_rapidity_difference_signal = [hist_rapidity_difference_signal.GetBinCenter(i) for i in range(1, hist_rapidity_difference_signal.GetNbinsX() + 1)]
y_vals_rapidity_difference_signal = [hist_rapidity_difference_signal.GetBinContent(i) for i in range(1, hist_rapidity_difference_signal.GetNbinsX() + 1)]

x_vals_rapidity_difference_background = [hist_rapidity_difference_background.GetBinCenter(i) for i in range(1, hist_rapidity_difference_background.GetNbinsX() + 1)]
y_vals_rapidity_difference_background = [hist_rapidity_difference_background.GetBinContent(i) for i in range(1, hist_rapidity_difference_background.GetNbinsX() + 1)]

# ✅ Convert ROOT histograms for missing transverse energy (MET)
x_vals_met_signal = [hist_met_signal.GetBinCenter(i) for i in range(1, hist_met_signal.GetNbinsX() + 1)]
y_vals_met_signal = [hist_met_signal.GetBinContent(i) for i in range(1, hist_met_signal.GetNbinsX() + 1)]

x_vals_met_background = [hist_met_background.GetBinCenter(i) for i in range(1, hist_met_background.GetNbinsX() + 1)]
y_vals_met_background = [hist_met_background.GetBinContent(i) for i in range(1, hist_met_background.GetNbinsX() + 1)]

# ✅ Convert ROOT histograms for subleading lepton pseudorapidity (η)
x_vals_subleading_lepton_eta_signal = [hist_subleading_lepton_eta_signal.GetBinCenter(i) for i in range(1, hist_subleading_lepton_eta_signal.GetNbinsX() + 1)]
y_vals_subleading_lepton_eta_signal = [hist_subleading_lepton_eta_signal.GetBinContent(i) for i in range(1, hist_subleading_lepton_eta_signal.GetNbinsX() + 1)]

x_vals_subleading_lepton_eta_background = [hist_subleading_lepton_eta_background.GetBinCenter(i) for i in range(1, hist_subleading_lepton_eta_background.GetNbinsX() + 1)]
y_vals_subleading_lepton_eta_background = [hist_subleading_lepton_eta_background.GetBinContent(i) for i in range(1, hist_subleading_lepton_eta_background.GetNbinsX() + 1)]





# ✅ Convert ROOT histograms for transverse mass of W boson (M_T(W))
x_vals_mT_W_signal = [hist_mT_W_signal.GetBinCenter(i) for i in range(1, hist_mT_W_signal.GetNbinsX() + 1)]
y_vals_mT_W_signal = [hist_mT_W_signal.GetBinContent(i) for i in range(1, hist_mT_W_signal.GetNbinsX() + 1)]

x_vals_mT_W_background = [hist_mT_W_background.GetBinCenter(i) for i in range(1, hist_mT_W_background.GetNbinsX() + 1)]
y_vals_mT_W_background = [hist_mT_W_background.GetBinContent(i) for i in range(1, hist_mT_W_background.GetNbinsX() + 1)]

# ✅ Convert ROOT histograms for reconstructed diboson system mass (M_WW)
x_vals_m_WW_signal = [hist_m_WW_signal.GetBinCenter(i) for i in range(1, hist_m_WW_signal.GetNbinsX() + 1)]
y_vals_m_WW_signal = [hist_m_WW_signal.GetBinContent(i) for i in range(1, hist_m_WW_signal.GetNbinsX() + 1)]

x_vals_m_WW_background = [hist_m_WW_background.GetBinCenter(i) for i in range(1, hist_m_WW_background.GetNbinsX() + 1)]
y_vals_m_WW_background = [hist_m_WW_background.GetBinContent(i) for i in range(1, hist_m_WW_background.GetNbinsX() + 1)]

# ✅ Convert ROOT histograms for diboson transverse momentum (pT_WW)
x_vals_pT_WW_signal = [hist_pT_WW_signal.GetBinCenter(i) for i in range(1, hist_pT_WW_signal.GetNbinsX() + 1)]
y_vals_pT_WW_signal = [hist_pT_WW_signal.GetBinContent(i) for i in range(1, hist_pT_WW_signal.GetNbinsX() + 1)]

x_vals_pT_WW_background = [hist_pT_WW_background.GetBinCenter(i) for i in range(1, hist_pT_WW_background.GetNbinsX() + 1)]
y_vals_pT_WW_background = [hist_pT_WW_background.GetBinContent(i) for i in range(1, hist_pT_WW_background.GetNbinsX() + 1)]

# ✅ Convert ROOT histograms for subleading lepton pseudorapidity (η)
x_vals_subleading_lepton_eta_signal = [hist_subleading_lepton_eta_signal.GetBinCenter(i) for i in range(1, hist_subleading_lepton_eta_signal.GetNbinsX() + 1)]
y_vals_subleading_lepton_eta_signal = [hist_subleading_lepton_eta_signal.GetBinContent(i) for i in range(1, hist_subleading_lepton_eta_signal.GetNbinsX() + 1)]

x_vals_subleading_lepton_eta_background = [hist_subleading_lepton_eta_background.GetBinCenter(i) for i in range(1, hist_subleading_lepton_eta_background.GetNbinsX() + 1)]
y_vals_subleading_lepton_eta_background = [hist_subleading_lepton_eta_background.GetBinContent(i) for i in range(1, hist_subleading_lepton_eta_background.GetNbinsX() + 1)]






#==============================================================================
#==============================================================================




# ✅ NEW: Plot the histograms for leading lepton pT
plt.figure(figsize=(10, 12))
plt.subplots_adjust(left=0.15, right=0.95, bottom=0.12, top=0.95)

plt.bar(x_vals_leading_lepton_signal, y_vals_leading_lepton_signal, width=hist_leading_lepton_signal.GetBinWidth(1), alpha=0.6, color="red", label="Signal ($W^+ W^-$) [$f_{M_2} / \Lambda^4$]")
plt.bar(x_vals_leading_lepton_background, y_vals_leading_lepton_background, width=hist_leading_lepton_background.GetBinWidth(1), alpha=0.6, color="blue", label="SM background ($W^+ W^-$)")
plt.xlabel(r"$p_T^{\mathrm{leading}~\ell} \ \mathrm{[GeV]}$")
plt.ylabel("Entries")
plt.title(r"Leading Lepton Transverse Momentum : $e^- p \to e^- W^+ W^- p \to e^- \ell^+ \nu_{\ell} \ell^- \bar{\nu}_{\ell} p$", fontsize=18)
plt.legend()
plt.grid()
plt.savefig("/home/hamzeh-khanpour/Documents/GitHub/LHeC_Fast_Simulation/Pythia8_Delphes_fully_leptonic/leading_lepton_pt_comparison.pdf", dpi=300)
plt.show()




# ✅ NEW: Plot the histograms for subleading lepton pT
plt.figure(figsize=(10, 12))
plt.subplots_adjust(left=0.15, right=0.95, bottom=0.12, top=0.95)

plt.bar(x_vals_subleading_lepton_signal, y_vals_subleading_lepton_signal, width=hist_subleading_lepton_signal.GetBinWidth(1), alpha=0.6, color="red", label="Signal ($W^+ W^-$) [$f_{M_2} / \Lambda^4$]")
plt.bar(x_vals_subleading_lepton_background, y_vals_subleading_lepton_background, width=hist_subleading_lepton_background.GetBinWidth(1), alpha=0.6, color="blue", label="SM background ($W^+ W^-$)")
plt.xlabel(r"$p_T^{\mathrm{subleading}~\ell} \ \mathrm{[GeV]}$")
plt.ylabel("Entries")
plt.title(r"Subleading Lepton Transverse Momentum : $e^- p \to e^- W^+ W^- p \to e^- \ell^+ \nu_{\ell} \ell^- \bar{\nu}_{\ell} p$", fontsize=18)
plt.legend()
plt.grid()
plt.savefig("/home/hamzeh-khanpour/Documents/GitHub/LHeC_Fast_Simulation/Pythia8_Delphes_fully_leptonic/subleading_lepton_pt_comparison.pdf", dpi=300)
plt.show()




# ✅ NEW: Plot the histograms for dilepton invariant mass
plt.figure(figsize=(10, 12))
plt.subplots_adjust(left=0.15, right=0.95, bottom=0.12, top=0.95)

plt.bar(x_vals_m_ll_signal, y_vals_m_ll_signal, width=hist_m_ll_signal.GetBinWidth(1), alpha=0.6, color="red", label="Signal ($W^+ W^-$) [$f_{M_2} / \Lambda^4$]")
plt.bar(x_vals_m_ll_background, y_vals_m_ll_background, width=hist_m_ll_background.GetBinWidth(1), alpha=0.6, color="blue", label="SM background ($W^+ W^-$)")
plt.xlabel(r"$M_{\ell^+\ell^-} \ \mathrm{[GeV]}$")
plt.ylabel("Entries")
plt.title(r"Dilepton Invariant Mass Distribution : $e^- p \to e^- W^+ W^- p \to e^- \ell^+ \nu_{\ell} \ell^- \bar{\nu}_{\ell} p$", fontsize=18)
plt.legend()
plt.grid()
plt.savefig("/home/hamzeh-khanpour/Documents/GitHub/LHeC_Fast_Simulation/Pythia8_Delphes_fully_leptonic/dilepton_mass_comparison.pdf", dpi=300)
plt.show()




# ✅ NEW: Plot the histograms for dilepton transverse momentum
plt.figure(figsize=(10, 12))
plt.subplots_adjust(left=0.15, right=0.95, bottom=0.12, top=0.95)

plt.bar(x_vals_pt_ll_signal, y_vals_pt_ll_signal, width=hist_pt_ll_signal.GetBinWidth(1), alpha=0.6, color="red", label="Signal ($W^+ W^-$) [$f_{M_2} / \Lambda^4$]")
plt.bar(x_vals_pt_ll_background, y_vals_pt_ll_background, width=hist_pt_ll_background.GetBinWidth(1), alpha=0.6, color="blue", label="SM background ($W^+ W^-$)")
plt.xlabel(r"$p_T^{\ell^+\ell^-} \ \mathrm{[GeV]}$")
plt.ylabel("Entries")
plt.title(r"Dilepton Transverse Momentum Distribution : $e^- p \to e^- W^+ W^- p \to e^- \ell^+ \nu_{\ell} \ell^- \bar{\nu}_{\ell} p$", fontsize=18)
plt.legend()
plt.grid()
plt.savefig("/home/hamzeh-khanpour/Documents/GitHub/LHeC_Fast_Simulation/Pythia8_Delphes_fully_leptonic/dilepton_pt_comparison.pdf", dpi=300)
plt.show()




# ✅ NEW: Plot the histograms for dilepton rapidity
plt.figure(figsize=(10, 12))
plt.subplots_adjust(left=0.15, right=0.95, bottom=0.12, top=0.95)

plt.bar(x_vals_rapidity_ll_signal, y_vals_rapidity_ll_signal, width=hist_rapidity_ll_signal.GetBinWidth(1), alpha=0.6, color="red", label="Signal ($W^+ W^-$) [$f_{M_2} / \Lambda^4$]")
plt.bar(x_vals_rapidity_ll_background, y_vals_rapidity_ll_background, width=hist_rapidity_ll_background.GetBinWidth(1), alpha=0.6, color="blue", label="SM background ($W^+ W^-$)")
plt.xlabel(r"$Y_{\ell^+\ell^-}$")
plt.ylabel("Entries")
plt.title(r"Dilepton System Rapidity : $e^- p \to e^- W^+ W^- p \to e^- \ell^+ \nu_{\ell} \ell^- \bar{\nu}_{\ell} p$", fontsize=18)
plt.legend()
plt.grid()
plt.savefig("/home/hamzeh-khanpour/Documents/GitHub/LHeC_Fast_Simulation/Pythia8_Delphes_fully_leptonic/dilepton_rapidity_comparison.pdf", dpi=300)
plt.show()




# ✅ NEW: Plot the histograms for leading lepton pseudorapidity
plt.figure(figsize=(10, 12))
plt.subplots_adjust(left=0.15, right=0.95, bottom=0.12, top=0.95)

plt.bar(x_vals_leading_lepton_eta_signal, y_vals_leading_lepton_eta_signal, width=hist_leading_lepton_eta_signal.GetBinWidth(1), alpha=0.6, color="red", label="Signal ($W^+ W^-$) [$f_{M_2} / \Lambda^4$]")
plt.bar(x_vals_leading_lepton_eta_background, y_vals_leading_lepton_eta_background, width=hist_leading_lepton_eta_background.GetBinWidth(1), alpha=0.6, color="blue", label="SM background ($W^+ W^-$)")
plt.xlabel(r"$\eta_{\mathrm{leading}~\ell}$")
plt.ylabel("Entries")
plt.title(r"Leading Lepton Pseudorapidity : $e^- p \to e^- W^+ W^- p \to e^- \ell^+ \nu_{\ell} \ell^- \bar{\nu}_{\ell} p$", fontsize=18)
plt.legend()
plt.grid()
plt.savefig("/home/hamzeh-khanpour/Documents/GitHub/LHeC_Fast_Simulation/Pythia8_Delphes_fully_leptonic/leading_lepton_eta_comparison.pdf", dpi=300)
plt.show()





# ✅ NEW: Plot the histograms for subleading lepton pseudorapidity
plt.figure(figsize=(10, 12))
plt.subplots_adjust(left=0.15, right=0.95, bottom=0.12, top=0.95)

plt.bar(x_vals_subleading_lepton_eta_signal, y_vals_subleading_lepton_eta_signal, width=hist_subleading_lepton_eta_signal.GetBinWidth(1), alpha=0.6, color="red", label="Signal ($W^+ W^-$) [$f_{M_2} / \Lambda^4$]")
plt.bar(x_vals_subleading_lepton_eta_background, y_vals_subleading_lepton_eta_background, width=hist_subleading_lepton_eta_background.GetBinWidth(1), alpha=0.6, color="blue", label="SM background ($W^+ W^-$)")
plt.xlabel(r"$\eta_{\mathrm{subleading}~\ell}$")
plt.ylabel("Entries")
plt.title(r"Subleading Lepton Pseudorapidity : $e^- p \to e^- W^+ W^- p \to e^- \ell^+ \nu_{\ell} \ell^- \bar{\nu}_{\ell} p$", fontsize=18)
plt.legend()
plt.grid()
plt.savefig("/home/hamzeh-khanpour/Documents/GitHub/LHeC_Fast_Simulation/Pythia8_Delphes_fully_leptonic/subleading_lepton_eta_comparison.pdf", dpi=300)
plt.show()





# ✅ NEW: Plot the histograms for lepton separation (ΔR)
plt.figure(figsize=(10, 12))
plt.subplots_adjust(left=0.15, right=0.95, bottom=0.12, top=0.95)

plt.bar(x_vals_lepton_separation_signal, y_vals_lepton_separation_signal, width=hist_lepton_separation_signal.GetBinWidth(1), alpha=0.6, color="red", label="Signal ($W^+ W^-$) [$f_{M_2} / \Lambda^4$]")
plt.bar(x_vals_lepton_separation_background, y_vals_lepton_separation_background, width=hist_lepton_separation_background.GetBinWidth(1), alpha=0.6, color="blue", label="SM background ($W^+ W^-$)")
plt.xlabel(r"$\Delta R(\ell_1, \ell_2)$")
plt.ylabel("Entries")
plt.title(r"Lepton Separation (ΔR) : $e^- p \to e^- W^+ W^- p \to e^- \ell^+ \nu_{\ell} \ell^- \bar{\nu}_{\ell} p$", fontsize=18)
plt.legend()
plt.grid()
plt.savefig("/home/hamzeh-khanpour/Documents/GitHub/LHeC_Fast_Simulation/Pythia8_Delphes_fully_leptonic/lepton_separation_comparison.pdf", dpi=300)
plt.show()





# ✅ NEW: Plot the histograms for azimuthal angle separation (Δφ)
plt.figure(figsize=(10, 12))
plt.subplots_adjust(left=0.15, right=0.95, bottom=0.12, top=0.95)

plt.bar(x_vals_azimuthal_angle_separation_signal, y_vals_azimuthal_angle_separation_signal, width=hist_azimuthal_angle_separation_signal.GetBinWidth(1), alpha=0.6, color="red", label="Signal ($W^+ W^-$) [$f_{M_2} / \Lambda^4$]")
plt.bar(x_vals_azimuthal_angle_separation_background, y_vals_azimuthal_angle_separation_background, width=hist_azimuthal_angle_separation_background.GetBinWidth(1), alpha=0.6, color="blue", label="SM background ($W^+ W^-$)")
plt.xlabel(r"$\Delta \phi (\ell_1, \ell_2)$")
plt.ylabel("Entries")
plt.title(r"Azimuthal Angle Separation (Δφ) : $e^- p \to e^- W^+ W^- p \to e^- \ell^+ \nu_{\ell} \ell^- \bar{\nu}_{\ell} p$", fontsize=18)
plt.legend()
plt.grid()
plt.savefig("/home/hamzeh-khanpour/Documents/GitHub/LHeC_Fast_Simulation/Pythia8_Delphes_fully_leptonic/azimuthal_angle_separation_comparison.pdf", dpi=300)
plt.show()





# ✅ NEW: Plot the histograms for rapidity difference (Δy)
plt.figure(figsize=(10, 12))
plt.subplots_adjust(left=0.15, right=0.95, bottom=0.12, top=0.95)

plt.bar(x_vals_rapidity_difference_signal, y_vals_rapidity_difference_signal, width=hist_rapidity_difference_signal.GetBinWidth(1), alpha=0.6, color="red", label="Signal ($W^+ W^-$) [$f_{M_2} / \Lambda^4$]")
plt.bar(x_vals_rapidity_difference_background, y_vals_rapidity_difference_background, width=hist_rapidity_difference_background.GetBinWidth(1), alpha=0.6, color="blue", label="SM background ($W^+ W^-$)")
plt.xlabel(r"$\Delta y (\ell_1, \ell_2)$")
plt.ylabel("Entries")
plt.title(r"Rapidity Difference (Δy) : $e^- p \to e^- W^+ W^- p \to e^- \ell^+ \nu_{\ell} \ell^- \bar{\nu}_{\ell} p$", fontsize=18)
plt.legend()
plt.grid()
plt.savefig("/home/hamzeh-khanpour/Documents/GitHub/LHeC_Fast_Simulation/Pythia8_Delphes_fully_leptonic/rapidity_difference_comparison.pdf", dpi=300)
plt.show()





# ✅ NEW: Plot the histograms for missing transverse energy (MET)
plt.figure(figsize=(10, 12))
plt.subplots_adjust(left=0.15, right=0.95, bottom=0.12, top=0.95)

plt.bar(x_vals_met_signal, y_vals_met_signal, width=hist_met_signal.GetBinWidth(1), alpha=0.6, color="red", label="Signal ($W^+ W^-$) [$f_{M_2} / \Lambda^4$]")
plt.bar(x_vals_met_background, y_vals_met_background, width=hist_met_background.GetBinWidth(1), alpha=0.6, color="blue", label="SM background ($W^+ W^-$)")
plt.xlabel(r"$E_T^{\mathrm{miss}}$ [GeV]")
plt.ylabel("Entries")
plt.title(r"Missing Transverse Energy (MET) : $e^- p \to e^- W^+ W^- p \to e^- \ell^+ \nu_{\ell} \ell^- \bar{\nu}_{\ell} p$", fontsize=18)
plt.legend()
plt.grid()
plt.savefig("/home/hamzeh-khanpour/Documents/GitHub/LHeC_Fast_Simulation/Pythia8_Delphes_fully_leptonic/met_comparison.pdf", dpi=300)
plt.show()






# ✅ NEW: Plot the histograms for transverse mass of W boson (M_T(W))
plt.figure(figsize=(10, 12))
plt.subplots_adjust(left=0.15, right=0.95, bottom=0.12, top=0.95)

plt.bar(x_vals_mT_W_signal, y_vals_mT_W_signal, width=hist_mT_W_signal.GetBinWidth(1), alpha=0.6, color="red", label="Signal ($W^+ W^-$) [$f_{M_2} / \Lambda^4$]")
plt.bar(x_vals_mT_W_background, y_vals_mT_W_background, width=hist_mT_W_background.GetBinWidth(1), alpha=0.6, color="blue", label="SM background ($W^+ W^-$)")
plt.xlabel(r"$M_T(W)$ [GeV]")
plt.ylabel("Entries")
plt.title(r"Transverse Mass of W Boson : $e^- p \to e^- W^+ W^- p \to e^- \ell^+ \nu_{\ell} \ell^- \bar{\nu}_{\ell} p$", fontsize=18)
plt.legend()
plt.grid()
plt.savefig("/home/hamzeh-khanpour/Documents/GitHub/LHeC_Fast_Simulation/Pythia8_Delphes_fully_leptonic/mT_W_comparison.pdf", dpi=300)
plt.show()





# ✅ NEW: Plot the histograms for reconstructed diboson system mass (M_WW)
plt.figure(figsize=(10, 12))
plt.subplots_adjust(left=0.15, right=0.95, bottom=0.12, top=0.95)

plt.bar(x_vals_m_WW_signal, y_vals_m_WW_signal, width=hist_m_WW_signal.GetBinWidth(1), alpha=0.6, color="red", label="Signal ($W^+ W^-$) [$f_{M_2} / \Lambda^4$]")
plt.bar(x_vals_m_WW_background, y_vals_m_WW_background, width=hist_m_WW_background.GetBinWidth(1), alpha=0.6, color="blue", label="SM background ($W^+ W^-$)")
plt.xlabel(r"$M_{WW}$ [GeV]")
plt.ylabel("Entries")
plt.title(r"Reconstructed Diboson System Mass : $e^- p \to e^- W^+ W^- p \to e^- \ell^+ \nu_{\ell} \ell^- \bar{\nu}_{\ell} p$", fontsize=18)
plt.legend()
plt.grid()
plt.savefig("/home/hamzeh-khanpour/Documents/GitHub/LHeC_Fast_Simulation/Pythia8_Delphes_fully_leptonic/m_WW_comparison.pdf", dpi=300)
plt.show()




# ✅ NEW: Plot the histograms for diboson transverse momentum (pT_WW)
plt.figure(figsize=(10, 12))
plt.subplots_adjust(left=0.15, right=0.95, bottom=0.12, top=0.95)

plt.bar(x_vals_pT_WW_signal, y_vals_pT_WW_signal, width=hist_pT_WW_signal.GetBinWidth(1), alpha=0.6, color="red", label="Signal ($W^+ W^-$) [$f_{M_2} / \Lambda^4$]")
plt.bar(x_vals_pT_WW_background, y_vals_pT_WW_background, width=hist_pT_WW_background.GetBinWidth(1), alpha=0.6, color="blue", label="SM background ($W^+ W^-$)")
plt.xlabel(r"$p_T^{WW}$ [GeV]")
plt.ylabel("Entries")
plt.title(r"Diboson Transverse Momentum : $e^- p \to e^- W^+ W^- p \to e^- \ell^+ \nu_{\ell} \ell^- \bar{\nu}_{\ell} p$", fontsize=18)
plt.legend()
plt.grid()
plt.savefig("/home/hamzeh-khanpour/Documents/GitHub/LHeC_Fast_Simulation/Pythia8_Delphes_fully_leptonic/pT_WW_comparison.pdf", dpi=300)
plt.show()















