

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
#    "$FM_{0} / \Lambda^4$": "aa_ww_semi_leptonic_NP_1_FM0_100.root",
#    "$FM_{1} / \Lambda^4$": "aa_ww_semi_leptonic_NP_1_FM1_100.root",
    "$FM_{2} / \Lambda^4$": "aa_ww_semi_leptonic_NP_1_FM2_100.root",
#    "$FM_{3} / \Lambda^4$": "aa_ww_semi_leptonic_NP_1_FM3_100.root",
}
aa_ww_background_file_path = "aa_ww_semi_leptonic_SM_decay.root"

aa_ttbar_background_file_path  = "aa_ttbar_inclusive_decay.root"
aa_tautau_background_file_path = "aa_tautau_leptonic_decay.root"

aa_mumu_background_file_path   = "aa_mumu.root"



inclusive_ttbar_background_file_path   = "inclusive_ttbar_inclusive_decay.root"
single_top_background_file_path        = "single_top_inclusive_decay.root"

w_production_background_file_path      = "w_production_inclusive_decay.root"
z_production_background_file_path      = "z_production_inclusive_decay.root"


wwj_production_background_file_path      = "diboson_wwj_inclusive_decay.root"
zzj_production_background_file_path      = "diboson_zzj_inclusive_decay.root"
wzj_production_background_file_path      = "diboson_wzj_inclusive_decay.root"


# Load Delphes library
ROOT.gSystem.Load("libDelphes")
ROOT.gInterpreter.Declare('#include "classes/DelphesClasses.h"')
ROOT.gInterpreter.Declare('#include "external/ExRootAnalysis/ExRootTreeReader.h"')


#=========================================================================
#=========================================================================


def process_file(
    file_path,
    hist_lepton_pt,
    hist_leading_jet_pt,
    hist_lepton_eta,
    hist_delta_r,
    hist_missing_et,
    hist_subleading_jet_eta,
    hist_leading_jet_eta,
    hist_jet_centrality,
    hist_delta_eta_jj,
    hist_m_w_leptonic,
    hist_m_w_hadronic,
):

    # Open the ROOT file
    chain = ROOT.TChain("Delphes")
    chain.Add(file_path)

    # Create ExRootTreeReader object
    treeReader = ROOT.ExRootTreeReader(chain)
    numberOfEntries = treeReader.GetEntries()

    # Counters for efficiency calculation
    # ✅ Initialize selection counters
    total_events = numberOfEntries

    selected_events_pre = 0   # ✅ Fix: Initialize before use
    selected_events_final = 0 # ✅ Fix: Initialize before use


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
        selected_events_pre += 1


        # Fill histogram with lepton pT
        hist_lepton_pt.Fill(leptons[0].Pt())
        hist_lepton_eta.Fill(leptons[0].Eta())   # ✅ Fill η histogram

        # Fill histogram with leading jet pT
        leading_jet = jets[0] if jets[0].Pt() > jets[1].Pt() else jets[1]
        hist_leading_jet_pt.Fill(leading_jet.Pt())


        # Fill optional histograms  Fill histogram with leading jet pT
        if hist_delta_r is not None:
            delta_r = leptons[0].DeltaR(leading_jet)
            hist_delta_r.Fill(delta_r)


        if hist_missing_et is not None and branchMissingET.GetEntries() > 0:
            missing_et = branchMissingET.At(0)
            hist_missing_et.Fill(missing_et.MET)


        if hist_subleading_jet_eta is not None and hist_leading_jet_eta is not None:
            delta_eta_jj = abs(jets[0].Eta() - jets[1].Eta())
            if delta_eta_jj > -100000:
                subleading_jet_eta = jets[1].Eta()
                hist_subleading_jet_eta.Fill(subleading_jet_eta)
                hist_leading_jet_eta.Fill(jets[0].Eta())


        if hist_jet_centrality is not None:
            jet_centrality = abs(jets[0].Eta() + jets[1].Eta()) / 2.0
            hist_jet_centrality.Fill(jet_centrality)


        if hist_delta_eta_jj is not None:
            delta_eta_jj = abs(jets[0].Eta() - jets[1].Eta())
            hist_delta_eta_jj.Fill(delta_eta_jj)


        # Calculate invariant mass of the hadronic W boson
        if hist_m_w_hadronic is not None:
            w_hadronic = jets[0] + jets[1]
            hist_m_w_hadronic.Fill(w_hadronic.M())



        # **✅ Corrected W → ℓν Reconstruction**
        if hist_m_w_leptonic is not None and branchMissingET.GetEntries() > 0:
            missing_et = branchMissingET.At(0)

            # **Construct MET components**
            px_nu = missing_et.MET * np.cos(missing_et.Phi)
            py_nu = missing_et.MET * np.sin(missing_et.Phi)


            # **W Boson Mass Constraint**
            MW = 80.385  # W boson mass in GeV
            A1 = 2 * (leptons[0].Px() * px_nu + leptons[0].Py() * py_nu) + MW**2
            B1 = 4 * leptons[0].E()**2 * (px_nu**2 + py_nu**2)
            C1 = 4 * (leptons[0].E()**2 - leptons[0].Pz()**2)
            D1 = -A1 * 4 * leptons[0].Pz()
            E1 = B1 - (A1**2) / 4  # ✅ Corrected denominator

            # **Solve for Pz using the quadratic equation**
            discriminant = D1**2 - 4 * C1 * E1
            if discriminant >= 0 and C1 != 0:
                Pz1 = (-D1 + np.sqrt(discriminant)) / (2 * C1)
                Pz2 = (-D1 - np.sqrt(discriminant)) / (2 * C1)
                Pz = Pz1 if abs(Pz1) < abs(Pz2) else Pz2  # Choose physical solution
            elif C1 != 0:
                Pz = -D1 / (2 * C1)  # Use fallback solution
            else:
                Pz = -E1 / D1 if D1 != 0 else 0  # Handle special case to avoid division by zero


            # **Construct final neutrino 4-vector with corrected energy**
            neutrino_vec = TLorentzVector()

            E_nu = np.sqrt(px_nu**2 + py_nu**2 + Pz**2)  # ✅ Corrected energy calculation
            neutrino_vec.SetPxPyPzE(px_nu, py_nu, Pz, E_nu)  # ✅ Corrected Pz


            # **Reconstruct the W boson correctly**
            w_leptonic = leptons[0] + neutrino_vec
            hist_m_w_leptonic.Fill(w_leptonic.M())


        # **✅ Final Event Selection (Corrected Indentation)**
        if leading_jet.Pt() < 1 or leptons[0].Pt() < 1:
            continue

        # Count selected events
        selected_events_final += 1



    # **✅ Selection Efficiency Calculations**
    efficiency_pre = selected_events_pre / total_events if total_events > 0 else 0

    efficiency_final = selected_events_final / total_events if total_events > 0 else 0



    # Return histograms and efficiency
    return (
        hist_lepton_pt,
        hist_leading_jet_pt,
        hist_lepton_eta,
        hist_delta_r,
        hist_missing_et,
        hist_subleading_jet_eta,
        hist_leading_jet_eta,
        hist_jet_centrality,
        hist_delta_eta_jj,
        hist_m_w_leptonic,
        hist_m_w_hadronic,
        efficiency_pre,
        efficiency_final,
    )






#=========================================================================
#=========================================================================



# Parameters for differential cross-section

#10
#signal_cross_sections = {
#    "$FM_{0} / \Lambda^4$": 0.0,   # pb
#    "$FM_{1} / \Lambda^4$": 0.0,   # pb
#    "$FM_{2} / \Lambda^4$": 0.0,   # pb
#    "$FM_{3} / \Lambda^4$": 0.0    # pb
#}



#100
signal_cross_sections = {
    "$FM_{0} / \Lambda^4$": 0.01490319,   # pb
    "$FM_{1} / \Lambda^4$": 0.01508150,   # pb
    "$FM_{2} / \Lambda^4$": 0.02142500,   # pb
    "$FM_{3} / \Lambda^4$": 0.01644609    # pb
}


aa_ww_background_cross_section = 0.0150743  # pb


aa_ttbar_background_cross_section  = 4.824851e-05 * 100.0  # pb  * 10^{+2}
aa_tautau_background_cross_section = 2.51510000  # pb
aa_mumu_background_cross_section   = 2.57270000  # pb


inclusive_ttbar_background_cross_section   = 0.0065764            # pb
single_top_background_cross_section        = 1.36209              # pb

w_production_background_cross_section      = 1.910288             # pb
z_production_background_cross_section      = 0.24064758729900002  # pb


wwj_production_background_cross_section      = 0.016080595320336195   # pb
zzj_production_background_cross_section      = 6.694889944457796e-05 * 100.0  # pb  * 10^{+2}
wzj_production_background_cross_section      = 0.0023785292894910495  # pb


num_bins = 50
pt_range_lepton = (0, 300)     # Range for lepton pT
pt_range_jet = (0, 300)        # Range for leading jet pT (adjusted for higher jet momenta)
eta_range = (-4, 6)          # Range for pseudorapidity
delta_r_range = (0, 6)        # Range for Delta R
met_range = (0, 300)           # Range for Missing Transverse Energy (MET)
centrality_range = (-3, 6)     # Range for centrality
exp_centrality_range = (-3, 6)  # Range for exponential centrality
jet_centrality_range = (-1, 6)  # Range for jet centrality
delta_eta_jj_range = (0, 6)    # Range for Delta Eta between jets


# Define ranges for invariant masses
m_w_hadronic_range = (1, 140)  # Range for the hadronic W boson mass
m_w_leptonic_range = (1, 140)  # Range for the leptonic W boson mass



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

# Bin width for W boson invariant masses
bin_width_m_w_hadronic = (m_w_hadronic_range[1] - m_w_hadronic_range[0]) / num_bins
bin_width_m_w_leptonic = (m_w_leptonic_range[1] - m_w_leptonic_range[0]) / num_bins




"""
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



def calculate_dsigma(histogram, total_cross_section, bin_width, verbose=False):
    """
    Compute differential cross section dσ/dX for a histogram.

    Args:
        histogram (ROOT.TH1): Input histogram with raw event counts.
        total_cross_section (float): Total cross section of the process in pb.
        bin_width (float): Bin width for normalization.
        verbose (bool): Print warning if histogram is empty.

    Returns:
        tuple: (list of bin edges, list of normalized dσ/dX values)
    """
    counts = [histogram.GetBinContent(i) for i in range(1, histogram.GetNbinsX() + 1)]
    bin_edges = [histogram.GetBinLowEdge(i) for i in range(1, histogram.GetNbinsX() + 2)]

    total_entries = sum(counts)

    if total_entries == 0:
        if verbose:
            print(f"⚠️ Warning: Histogram '{histogram.GetName()}' is empty. Skipping dσ/dX calculation.")
        dsigma = [0.0] * len(counts)
    else:
        norm = total_cross_section / total_entries / bin_width
        dsigma = [count * norm for count in counts]

    return bin_edges[:-1], dsigma






#=========================================================================
#=========================================================================



# Create histograms for lepton and leading jet pT
# Dictionary to store histograms for each signal
signal_histograms = {}

for signal_name in signal_files:
    signal_histograms[signal_name] = {
        "hist_lepton_pt": ROOT.TH1F(f"hist_lepton_{signal_name}", f"Lepton pT Distribution ({signal_name}); p_{{T}} [GeV]; Entries", num_bins, *pt_range_lepton),
        "hist_leading_jet_pt": ROOT.TH1F(f"hist_leading_jet_{signal_name}", f"Leading Jet pT Distribution ({signal_name}); p_{{T}} [GeV]; Entries", num_bins, *pt_range_jet),
        "hist_lepton_eta": ROOT.TH1F(f"hist_lepton_eta_{signal_name}", f"Lepton Eta Distribution ({signal_name}); #eta; Entries", num_bins, *eta_range),
        "hist_delta_r": ROOT.TH1F(f"hist_delta_r_{signal_name}", f"Delta R Distribution ({signal_name}); ΔR; Entries", num_bins, *delta_r_range),
        "hist_missing_et": ROOT.TH1F(f"hist_missing_et_{signal_name}", f"Missing ET Distribution ({signal_name}); MET [GeV]; Entries", num_bins, *met_range),
        "hist_subleading_jet_eta": ROOT.TH1F(f"hist_subleading_jet_eta_{signal_name}", f"Centrality Distribution ({signal_name}); Centrality; Entries", num_bins, *centrality_range),
        "hist_leading_jet_eta": ROOT.TH1F(f"hist_leading_jet_eta_{signal_name}", f"Exponential Centrality Distribution ({signal_name}); exp(Centrality); Entries", num_bins, *exp_centrality_range),
        "hist_jet_centrality": ROOT.TH1F(f"hist_jet_centrality_{signal_name}", f"Jet Centrality Distribution ({signal_name}); Jet Centrality; Entries", num_bins, *jet_centrality_range),
        "hist_delta_eta_jj": ROOT.TH1F(f"hist_delta_eta_jj_{signal_name}", f"Delta Eta Between Jets ({signal_name}); Δηjj; Entries", num_bins, *delta_eta_jj_range),
        "hist_m_w_leptonic": ROOT.TH1F(f"hist_m_w_leptonic_{signal_name}", f"Leptonic W Boson Mass ({signal_name}); m_{{W}}^{{leptonic}} [GeV]; Entries", num_bins, *m_w_leptonic_range),
        "hist_m_w_hadronic": ROOT.TH1F(f"hist_m_w_hadronic_{signal_name}", f"Hadronic W Boson Mass ({signal_name}); m_{{W}}^{{hadronic}} [GeV]; Entries", num_bins, *m_w_hadronic_range)
    }






# === AA_WW Background ===
hist_lepton_aa_ww = ROOT.TH1F("hist_lepton_aa_ww", "Lepton pT (aa_ww); p_{T} [GeV]; Entries", num_bins, *pt_range_lepton)
hist_leading_jet_aa_ww = ROOT.TH1F("hist_leading_jet_aa_ww", "Leading Jet pT (aa_ww); p_{T} [GeV]; Entries", num_bins, *pt_range_jet)
hist_lepton_eta_aa_ww = ROOT.TH1F("hist_lepton_eta_aa_ww", "Lepton Eta (aa_ww); #eta; Entries", num_bins, *eta_range)
hist_delta_r_aa_ww = ROOT.TH1F("hist_delta_r_aa_ww", "Delta R (aa_ww); ΔR; Entries", num_bins, *delta_r_range)
hist_missing_et_aa_ww = ROOT.TH1F("hist_missing_et_aa_ww", "MET (aa_ww); MET [GeV]; Entries", num_bins, *met_range)
hist_subleading_jet_eta_aa_ww = ROOT.TH1F("hist_subleading_jet_eta_aa_ww", "Centrality (aa_ww); Centrality; Entries", num_bins, *centrality_range)
hist_leading_jet_eta_aa_ww = ROOT.TH1F("hist_leading_jet_eta_aa_ww", "Exp Centrality (aa_ww); exp(Centrality); Entries", num_bins, *exp_centrality_range)
hist_jet_centrality_aa_ww = ROOT.TH1F("hist_jet_centrality_aa_ww", "Jet Centrality (aa_ww); Jet Centrality; Entries", num_bins, *jet_centrality_range)
hist_delta_eta_jj_aa_ww = ROOT.TH1F("hist_delta_eta_jj_aa_ww", "Delta Eta jj (aa_ww); Δηjj; Entries", num_bins, *delta_eta_jj_range)
hist_m_w_hadronic_aa_ww = ROOT.TH1F("hist_m_w_hadronic_aa_ww", "Hadronic W Mass (aa_ww); m_{W}^{hadronic} [GeV]; Entries", num_bins, *m_w_hadronic_range)
hist_m_w_leptonic_aa_ww = ROOT.TH1F("hist_m_w_leptonic_aa_ww", "Leptonic W Mass (aa_ww); m_{W}^{leptonic} [GeV]; Entries", num_bins, *m_w_leptonic_range)

# === AA_TTBAR ===
hist_lepton_aa_ttbar = ROOT.TH1F("hist_lepton_aa_ttbar", "Lepton pT (aa_ttbar); p_{T} [GeV]; Entries", num_bins, *pt_range_lepton)
hist_leading_jet_aa_ttbar = ROOT.TH1F("hist_leading_jet_aa_ttbar", "Leading Jet pT (aa_ttbar); p_{T} [GeV]; Entries", num_bins, *pt_range_jet)
hist_lepton_eta_aa_ttbar = ROOT.TH1F("hist_lepton_eta_aa_ttbar", "Lepton Eta (aa_ttbar); #eta; Entries", num_bins, *eta_range)
hist_delta_r_aa_ttbar = ROOT.TH1F("hist_delta_r_aa_ttbar", "Delta R (aa_ttbar); ΔR; Entries", num_bins, *delta_r_range)
hist_missing_et_aa_ttbar = ROOT.TH1F("hist_missing_et_aa_ttbar", "MET (aa_ttbar); MET [GeV]; Entries", num_bins, *met_range)
hist_subleading_jet_eta_aa_ttbar = ROOT.TH1F("hist_subleading_jet_eta_aa_ttbar", "Centrality (aa_ttbar); Centrality; Entries", num_bins, *centrality_range)
hist_leading_jet_eta_aa_ttbar = ROOT.TH1F("hist_leading_jet_eta_aa_ttbar", "Exp Centrality (aa_ttbar); exp(Centrality); Entries", num_bins, *exp_centrality_range)
hist_jet_centrality_aa_ttbar = ROOT.TH1F("hist_jet_centrality_aa_ttbar", "Jet Centrality (aa_ttbar); Jet Centrality; Entries", num_bins, *jet_centrality_range)
hist_delta_eta_jj_aa_ttbar = ROOT.TH1F("hist_delta_eta_jj_aa_ttbar", "Delta Eta jj (aa_ttbar); Δηjj; Entries", num_bins, *delta_eta_jj_range)
hist_m_w_hadronic_aa_ttbar = ROOT.TH1F("hist_m_w_hadronic_aa_ttbar", "Hadronic W Mass (aa_ttbar); m_{W}^{hadronic} [GeV]; Entries", num_bins, *m_w_hadronic_range)
hist_m_w_leptonic_aa_ttbar = ROOT.TH1F("hist_m_w_leptonic_aa_ttbar", "Leptonic W Mass (aa_ttbar); m_{W}^{leptonic} [GeV]; Entries", num_bins, *m_w_leptonic_range)

# === AA_TAUTAU ===
hist_lepton_aa_tautau = ROOT.TH1F("hist_lepton_aa_tautau", "Lepton pT (aa_tautau); p_{T} [GeV]; Entries", num_bins, *pt_range_lepton)
hist_leading_jet_aa_tautau = ROOT.TH1F("hist_leading_jet_aa_tautau", "Leading Jet pT (aa_tautau); p_{T} [GeV]; Entries", num_bins, *pt_range_jet)
hist_lepton_eta_aa_tautau = ROOT.TH1F("hist_lepton_eta_aa_tautau", "Lepton Eta (aa_tautau); #eta; Entries", num_bins, *eta_range)
hist_delta_r_aa_tautau = ROOT.TH1F("hist_delta_r_aa_tautau", "Delta R (aa_tautau); ΔR; Entries", num_bins, *delta_r_range)
hist_missing_et_aa_tautau = ROOT.TH1F("hist_missing_et_aa_tautau", "MET (aa_tautau); MET [GeV]; Entries", num_bins, *met_range)
hist_subleading_jet_eta_aa_tautau = ROOT.TH1F("hist_subleading_jet_eta_aa_tautau", "Centrality (aa_tautau); Centrality; Entries", num_bins, *centrality_range)
hist_leading_jet_eta_aa_tautau = ROOT.TH1F("hist_leading_jet_eta_aa_tautau", "Exp Centrality (aa_tautau); exp(Centrality); Entries", num_bins, *exp_centrality_range)
hist_jet_centrality_aa_tautau = ROOT.TH1F("hist_jet_centrality_aa_tautau", "Jet Centrality (aa_tautau); Jet Centrality; Entries", num_bins, *jet_centrality_range)
hist_delta_eta_jj_aa_tautau = ROOT.TH1F("hist_delta_eta_jj_aa_tautau", "Delta Eta jj (aa_tautau); Δηjj; Entries", num_bins, *delta_eta_jj_range)
hist_m_w_hadronic_aa_tautau = ROOT.TH1F("hist_m_w_hadronic_aa_tautau", "Hadronic W Mass (aa_tautau); m_{W}^{hadronic} [GeV]; Entries", num_bins, *m_w_hadronic_range)
hist_m_w_leptonic_aa_tautau = ROOT.TH1F("hist_m_w_leptonic_aa_tautau", "Leptonic W Mass (aa_tautau); m_{W}^{leptonic} [GeV]; Entries", num_bins, *m_w_leptonic_range)

# === AA_MUMU ===
hist_lepton_aa_mumu = ROOT.TH1F("hist_lepton_aa_mumu", "Lepton pT (aa_mumu); p_{T} [GeV]; Entries", num_bins, *pt_range_lepton)
hist_leading_jet_aa_mumu = ROOT.TH1F("hist_leading_jet_aa_mumu", "Leading Jet pT (aa_mumu); p_{T} [GeV]; Entries", num_bins, *pt_range_jet)
hist_lepton_eta_aa_mumu = ROOT.TH1F("hist_lepton_eta_aa_mumu", "Lepton Eta (aa_mumu); #eta; Entries", num_bins, *eta_range)
hist_delta_r_aa_mumu = ROOT.TH1F("hist_delta_r_aa_mumu", "Delta R (aa_mumu); ΔR; Entries", num_bins, *delta_r_range)
hist_missing_et_aa_mumu = ROOT.TH1F("hist_missing_et_aa_mumu", "MET (aa_mumu); MET [GeV]; Entries", num_bins, *met_range)
hist_subleading_jet_eta_aa_mumu = ROOT.TH1F("hist_subleading_jet_eta_aa_mumu", "Centrality (aa_mumu); Centrality; Entries", num_bins, *centrality_range)
hist_leading_jet_eta_aa_mumu = ROOT.TH1F("hist_leading_jet_eta_aa_mumu", "Exp Centrality (aa_mumu); exp(Centrality); Entries", num_bins, *exp_centrality_range)
hist_jet_centrality_aa_mumu = ROOT.TH1F("hist_jet_centrality_aa_mumu", "Jet Centrality (aa_mumu); Jet Centrality; Entries", num_bins, *jet_centrality_range)
hist_delta_eta_jj_aa_mumu = ROOT.TH1F("hist_delta_eta_jj_aa_mumu", "Delta Eta jj (aa_mumu); Δηjj; Entries", num_bins, *delta_eta_jj_range)
hist_m_w_hadronic_aa_mumu = ROOT.TH1F("hist_m_w_hadronic_aa_mumu", "Hadronic W Mass (aa_mumu); m_{W}^{hadronic} [GeV]; Entries", num_bins, *m_w_hadronic_range)
hist_m_w_leptonic_aa_mumu = ROOT.TH1F("hist_m_w_leptonic_aa_mumu", "Leptonic W Mass (aa_mumu); m_{W}^{leptonic} [GeV]; Entries", num_bins, *m_w_leptonic_range)

# === Inclusive TTBAR ===
hist_lepton_inclusive_ttbar = ROOT.TH1F("hist_lepton_inclusive_ttbar", "Lepton pT (inclusive_ttbar); p_{T} [GeV]; Entries", num_bins, *pt_range_lepton)
hist_leading_jet_inclusive_ttbar = ROOT.TH1F("hist_leading_jet_inclusive_ttbar", "Leading Jet pT (inclusive_ttbar); p_{T} [GeV]; Entries", num_bins, *pt_range_jet)
hist_lepton_eta_inclusive_ttbar = ROOT.TH1F("hist_lepton_eta_inclusive_ttbar", "Lepton Eta (inclusive_ttbar); #eta; Entries", num_bins, *eta_range)
hist_delta_r_inclusive_ttbar = ROOT.TH1F("hist_delta_r_inclusive_ttbar", "Delta R (inclusive_ttbar); ΔR; Entries", num_bins, *delta_r_range)
hist_missing_et_inclusive_ttbar = ROOT.TH1F("hist_missing_et_inclusive_ttbar", "MET (inclusive_ttbar); MET [GeV]; Entries", num_bins, *met_range)
hist_subleading_jet_eta_inclusive_ttbar = ROOT.TH1F("hist_subleading_jet_eta_inclusive_ttbar", "Centrality (inclusive_ttbar); Centrality; Entries", num_bins, *centrality_range)
hist_leading_jet_eta_inclusive_ttbar = ROOT.TH1F("hist_leading_jet_eta_inclusive_ttbar", "Exp Centrality (inclusive_ttbar); exp(Centrality); Entries", num_bins, *exp_centrality_range)
hist_jet_centrality_inclusive_ttbar = ROOT.TH1F("hist_jet_centrality_inclusive_ttbar", "Jet Centrality (inclusive_ttbar); Jet Centrality; Entries", num_bins, *jet_centrality_range)
hist_delta_eta_jj_inclusive_ttbar = ROOT.TH1F("hist_delta_eta_jj_inclusive_ttbar", "Delta Eta jj (inclusive_ttbar); Δηjj; Entries", num_bins, *delta_eta_jj_range)
hist_m_w_hadronic_inclusive_ttbar = ROOT.TH1F("hist_m_w_hadronic_inclusive_ttbar", "Hadronic W Mass (inclusive_ttbar); m_{W}^{hadronic} [GeV]; Entries", num_bins, *m_w_hadronic_range)
hist_m_w_leptonic_inclusive_ttbar = ROOT.TH1F("hist_m_w_leptonic_inclusive_ttbar", "Leptonic W Mass (inclusive_ttbar); m_{W}^{leptonic} [GeV]; Entries", num_bins, *m_w_leptonic_range)

# === SINGLE TOP ===
hist_lepton_single_top = ROOT.TH1F("hist_lepton_single_top", "Lepton pT (single_top); p_{T} [GeV]; Entries", num_bins, *pt_range_lepton)
hist_leading_jet_single_top = ROOT.TH1F("hist_leading_jet_single_top", "Leading Jet pT (single_top); p_{T} [GeV]; Entries", num_bins, *pt_range_jet)
hist_lepton_eta_single_top = ROOT.TH1F("hist_lepton_eta_single_top", "Lepton Eta (single_top); #eta; Entries", num_bins, *eta_range)
hist_delta_r_single_top = ROOT.TH1F("hist_delta_r_single_top", "Delta R (single_top); ΔR; Entries", num_bins, *delta_r_range)
hist_missing_et_single_top = ROOT.TH1F("hist_missing_et_single_top", "MET (single_top); MET [GeV]; Entries", num_bins, *met_range)
hist_subleading_jet_eta_single_top = ROOT.TH1F("hist_subleading_jet_eta_single_top", "Centrality (single_top); Centrality; Entries", num_bins, *centrality_range)
hist_leading_jet_eta_single_top = ROOT.TH1F("hist_leading_jet_eta_single_top", "Exp Centrality (single_top); exp(Centrality); Entries", num_bins, *exp_centrality_range)
hist_jet_centrality_single_top = ROOT.TH1F("hist_jet_centrality_single_top", "Jet Centrality (single_top); Jet Centrality; Entries", num_bins, *jet_centrality_range)
hist_delta_eta_jj_single_top = ROOT.TH1F("hist_delta_eta_jj_single_top", "Delta Eta jj (single_top); Δηjj; Entries", num_bins, *delta_eta_jj_range)
hist_m_w_hadronic_single_top = ROOT.TH1F("hist_m_w_hadronic_single_top", "Hadronic W Mass (single_top); m_{W}^{hadronic} [GeV]; Entries", num_bins, *m_w_hadronic_range)
hist_m_w_leptonic_single_top = ROOT.TH1F("hist_m_w_leptonic_single_top", "Leptonic W Mass (single_top); m_{W}^{leptonic} [GeV]; Entries", num_bins, *m_w_leptonic_range)

# === W PRODUCTION ===
hist_lepton_w_production = ROOT.TH1F("hist_lepton_w_production", "Lepton pT (w_production); p_{T} [GeV]; Entries", num_bins, *pt_range_lepton)
hist_leading_jet_w_production = ROOT.TH1F("hist_leading_jet_w_production", "Leading Jet pT (w_production); p_{T} [GeV]; Entries", num_bins, *pt_range_jet)
hist_lepton_eta_w_production = ROOT.TH1F("hist_lepton_eta_w_production", "Lepton Eta (w_production); #eta; Entries", num_bins, *eta_range)
hist_delta_r_w_production = ROOT.TH1F("hist_delta_r_w_production", "Delta R (w_production); ΔR; Entries", num_bins, *delta_r_range)
hist_missing_et_w_production = ROOT.TH1F("hist_missing_et_w_production", "MET (w_production); MET [GeV]; Entries", num_bins, *met_range)
hist_subleading_jet_eta_w_production = ROOT.TH1F("hist_subleading_jet_eta_w_production", "Centrality (w_production); Centrality; Entries", num_bins, *centrality_range)
hist_leading_jet_eta_w_production = ROOT.TH1F("hist_leading_jet_eta_w_production", "Exp Centrality (w_production); exp(Centrality); Entries", num_bins, *exp_centrality_range)
hist_jet_centrality_w_production = ROOT.TH1F("hist_jet_centrality_w_production", "Jet Centrality (w_production); Jet Centrality; Entries", num_bins, *jet_centrality_range)
hist_delta_eta_jj_w_production = ROOT.TH1F("hist_delta_eta_jj_w_production", "Delta Eta jj (w_production); Δηjj; Entries", num_bins, *delta_eta_jj_range)
hist_m_w_hadronic_w_production = ROOT.TH1F("hist_m_w_hadronic_w_production", "Hadronic W Mass (w_production); m_{W}^{hadronic} [GeV]; Entries", num_bins, *m_w_hadronic_range)
hist_m_w_leptonic_w_production = ROOT.TH1F("hist_m_w_leptonic_w_production", "Leptonic W Mass (w_production); m_{W}^{leptonic} [GeV]; Entries", num_bins, *m_w_leptonic_range)

# === Z PRODUCTION ===
hist_lepton_z_production = ROOT.TH1F("hist_lepton_z_production", "Lepton pT (z_production); p_{T} [GeV]; Entries", num_bins, *pt_range_lepton)
hist_leading_jet_z_production = ROOT.TH1F("hist_leading_jet_z_production", "Leading Jet pT (z_production); p_{T} [GeV]; Entries", num_bins, *pt_range_jet)
hist_lepton_eta_z_production = ROOT.TH1F("hist_lepton_eta_z_production", "Lepton Eta (z_production); #eta; Entries", num_bins, *eta_range)
hist_delta_r_z_production = ROOT.TH1F("hist_delta_r_z_production", "Delta R (z_production); ΔR; Entries", num_bins, *delta_r_range)
hist_missing_et_z_production = ROOT.TH1F("hist_missing_et_z_production", "MET (z_production); MET [GeV]; Entries", num_bins, *met_range)
hist_subleading_jet_eta_z_production = ROOT.TH1F("hist_subleading_jet_eta_z_production", "Centrality (z_production); Centrality; Entries", num_bins, *centrality_range)
hist_leading_jet_eta_z_production = ROOT.TH1F("hist_leading_jet_eta_z_production", "Exp Centrality (z_production); exp(Centrality); Entries", num_bins, *exp_centrality_range)
hist_jet_centrality_z_production = ROOT.TH1F("hist_jet_centrality_z_production", "Jet Centrality (z_production); Jet Centrality; Entries", num_bins, *jet_centrality_range)
hist_delta_eta_jj_z_production = ROOT.TH1F("hist_delta_eta_jj_z_production", "Delta Eta jj (z_production); Δηjj; Entries", num_bins, *delta_eta_jj_range)
hist_m_w_hadronic_z_production = ROOT.TH1F("hist_m_w_hadronic_z_production", "Hadronic W Mass (z_production); m_{W}^{hadronic} [GeV]; Entries", num_bins, *m_w_hadronic_range)
hist_m_w_leptonic_z_production = ROOT.TH1F("hist_m_w_leptonic_z_production", "Leptonic W Mass (z_production); m_{W}^{leptonic} [GeV]; Entries", num_bins, *m_w_leptonic_range)


# === WWJ PRODUCTION ===
hist_lepton_wwj_production = ROOT.TH1F("hist_lepton_wwj_production", "Lepton pT (wwj_production); p_{T} [GeV]; Entries", num_bins, *pt_range_lepton)
hist_leading_jet_wwj_production = ROOT.TH1F("hist_leading_jet_wwj_production", "Leading Jet pT (wwj_production); p_{T} [GeV]; Entries", num_bins, *pt_range_jet)
hist_lepton_eta_wwj_production = ROOT.TH1F("hist_lepton_eta_wwj_production", "Lepton Eta (wwj_production); #eta; Entries", num_bins, *eta_range)
hist_delta_r_wwj_production = ROOT.TH1F("hist_delta_r_wwj_production", "Delta R (wwj_production); ΔR; Entries", num_bins, *delta_r_range)
hist_missing_et_wwj_production = ROOT.TH1F("hist_missing_et_wwj_production", "MET (wwj_production); MET [GeV]; Entries", num_bins, *met_range)
hist_subleading_jet_eta_wwj_production = ROOT.TH1F("hist_subleading_jet_eta_wwj_production", "Centrality (wwj_production); Centrality; Entries", num_bins, *centrality_range)
hist_leading_jet_eta_wwj_production = ROOT.TH1F("hist_leading_jet_eta_wwj_production", "Exp Centrality (wwj_production); exp(Centrality); Entries", num_bins, *exp_centrality_range)
hist_jet_centrality_wwj_production = ROOT.TH1F("hist_jet_centrality_wwj_production", "Jet Centrality (wwj_production); Jet Centrality; Entries", num_bins, *jet_centrality_range)
hist_delta_eta_jj_wwj_production = ROOT.TH1F("hist_delta_eta_jj_wwj_production", "Delta Eta jj (wwj_production); Δηjj; Entries", num_bins, *delta_eta_jj_range)
hist_m_w_hadronic_wwj_production = ROOT.TH1F("hist_m_w_hadronic_wwj_production", "Hadronic W Mass (wwj_production); m_{W}^{hadronic} [GeV]; Entries", num_bins, *m_w_hadronic_range)
hist_m_w_leptonic_wwj_production = ROOT.TH1F("hist_m_w_leptonic_wwj_production", "Leptonic W Mass (wwj_production); m_{W}^{leptonic} [GeV]; Entries", num_bins, *m_w_leptonic_range)


# === ZZJ PRODUCTION ===
hist_lepton_zzj_production = ROOT.TH1F("hist_lepton_zzj_production", "Lepton pT (zzj_production); p_{T} [GeV]; Entries", num_bins, *pt_range_lepton)
hist_leading_jet_zzj_production = ROOT.TH1F("hist_leading_jet_zzj_production", "Leading Jet pT (zzj_production); p_{T} [GeV]; Entries", num_bins, *pt_range_jet)
hist_lepton_eta_zzj_production = ROOT.TH1F("hist_lepton_eta_zzj_production", "Lepton Eta (zzj_production); #eta; Entries", num_bins, *eta_range)
hist_delta_r_zzj_production = ROOT.TH1F("hist_delta_r_zzj_production", "Delta R (zzj_production); ΔR; Entries", num_bins, *delta_r_range)
hist_missing_et_zzj_production = ROOT.TH1F("hist_missing_et_zzj_production", "MET (zzj_production); MET [GeV]; Entries", num_bins, *met_range)
hist_subleading_jet_eta_zzj_production = ROOT.TH1F("hist_subleading_jet_eta_zzj_production", "Centrality (zzj_production); Centrality; Entries", num_bins, *centrality_range)
hist_leading_jet_eta_zzj_production = ROOT.TH1F("hist_leading_jet_eta_zzj_production", "Exp Centrality (zzj_production); exp(Centrality); Entries", num_bins, *exp_centrality_range)
hist_jet_centrality_zzj_production = ROOT.TH1F("hist_jet_centrality_zzj_production", "Jet Centrality (zzj_production); Jet Centrality; Entries", num_bins, *jet_centrality_range)
hist_delta_eta_jj_zzj_production = ROOT.TH1F("hist_delta_eta_jj_zzj_production", "Delta Eta jj (zzj_production); Δηjj; Entries", num_bins, *delta_eta_jj_range)
hist_m_w_hadronic_zzj_production = ROOT.TH1F("hist_m_w_hadronic_zzj_production", "Hadronic W Mass (zzj_production); m_{W}^{hadronic} [GeV]; Entries", num_bins, *m_w_hadronic_range)
hist_m_w_leptonic_zzj_production = ROOT.TH1F("hist_m_w_leptonic_zzj_production", "Leptonic W Mass (zzj_production); m_{W}^{leptonic} [GeV]; Entries", num_bins, *m_w_leptonic_range)


# === WZJ PRODUCTION ===
hist_lepton_wzj_production = ROOT.TH1F("hist_lepton_wzj_production", "Lepton pT (wzj_production); p_{T} [GeV]; Entries", num_bins, *pt_range_lepton)
hist_leading_jet_wzj_production = ROOT.TH1F("hist_leading_jet_wzj_production", "Leading Jet pT (wzj_production); p_{T} [GeV]; Entries", num_bins, *pt_range_jet)
hist_lepton_eta_wzj_production = ROOT.TH1F("hist_lepton_eta_wzj_production", "Lepton Eta (wzj_production); #eta; Entries", num_bins, *eta_range)
hist_delta_r_wzj_production = ROOT.TH1F("hist_delta_r_wzj_production", "Delta R (wzj_production); ΔR; Entries", num_bins, *delta_r_range)
hist_missing_et_wzj_production = ROOT.TH1F("hist_missing_et_wzj_production", "MET (wzj_production); MET [GeV]; Entries", num_bins, *met_range)
hist_subleading_jet_eta_wzj_production = ROOT.TH1F("hist_subleading_jet_eta_wzj_production", "Centrality (wzj_production); Centrality; Entries", num_bins, *centrality_range)
hist_leading_jet_eta_wzj_production = ROOT.TH1F("hist_leading_jet_eta_wzj_production", "Exp Centrality (wzj_production); exp(Centrality); Entries", num_bins, *exp_centrality_range)
hist_jet_centrality_wzj_production = ROOT.TH1F("hist_jet_centrality_wzj_production", "Jet Centrality (wzj_production); Jet Centrality; Entries", num_bins, *jet_centrality_range)
hist_delta_eta_jj_wzj_production = ROOT.TH1F("hist_delta_eta_jj_wzj_production", "Delta Eta jj (wzj_production); Δηjj; Entries", num_bins, *delta_eta_jj_range)
hist_m_w_hadronic_wzj_production = ROOT.TH1F("hist_m_w_hadronic_wzj_production", "Hadronic W Mass (wzj_production); m_{W}^{hadronic} [GeV]; Entries", num_bins, *m_w_hadronic_range)
hist_m_w_leptonic_wzj_production = ROOT.TH1F("hist_m_w_leptonic_wzj_production", "Leptonic W Mass (wzj_production); m_{W}^{leptonic} [GeV]; Entries", num_bins, *m_w_leptonic_range)










# Dictionary to store efficiencies for each signal
signal_efficiencies = {}

# Process all signal files dynamically
for signal_name, file_path in signal_files.items():
    print(f"Processing signal: {signal_name}")

    # Get the corresponding histograms for this signal
    histograms = signal_histograms[signal_name]

    # Process the signal file
    (histograms["hist_lepton_pt"], histograms["hist_leading_jet_pt"], histograms["hist_lepton_eta"], histograms["hist_delta_r"],
     histograms["hist_missing_et"], histograms["hist_subleading_jet_eta"], histograms["hist_leading_jet_eta"], histograms["hist_jet_centrality"],
     histograms["hist_delta_eta_jj"], histograms["hist_m_w_leptonic"], histograms["hist_m_w_hadronic"],
     efficiency_pre, efficiency_final) = process_file(
        file_path, histograms["hist_lepton_pt"], histograms["hist_leading_jet_pt"], histograms["hist_lepton_eta"],
        histograms["hist_delta_r"], histograms["hist_missing_et"], histograms["hist_subleading_jet_eta"], histograms["hist_leading_jet_eta"],
        histograms["hist_jet_centrality"], histograms["hist_delta_eta_jj"], histograms["hist_m_w_leptonic"], histograms["hist_m_w_hadronic"]
    )

    # Store efficiencies for this signal
    signal_efficiencies[signal_name] = {
        "efficiency_pre": efficiency_pre,
        "efficiency_final": efficiency_final
    }






# === Process the AA_WW background file ===
(hist_lepton_aa_ww, hist_leading_jet_aa_ww, hist_lepton_eta_aa_ww, hist_delta_r_aa_ww,
 hist_missing_et_aa_ww, hist_subleading_jet_eta_aa_ww, hist_leading_jet_eta_aa_ww, hist_jet_centrality_aa_ww,
 hist_delta_eta_jj_aa_ww, hist_m_w_leptonic_aa_ww, hist_m_w_hadronic_aa_ww,
 background_efficiency_pre_aa_ww, background_efficiency_final_aa_ww) = process_file(
    aa_ww_background_file_path, hist_lepton_aa_ww, hist_leading_jet_aa_ww, hist_lepton_eta_aa_ww,
    hist_delta_r_aa_ww, hist_missing_et_aa_ww, hist_subleading_jet_eta_aa_ww, hist_leading_jet_eta_aa_ww,
    hist_jet_centrality_aa_ww, hist_delta_eta_jj_aa_ww, hist_m_w_leptonic_aa_ww, hist_m_w_hadronic_aa_ww
)

# === Process the AA_TTBAR background file ===
(hist_lepton_aa_ttbar, hist_leading_jet_aa_ttbar, hist_lepton_eta_aa_ttbar, hist_delta_r_aa_ttbar,
 hist_missing_et_aa_ttbar, hist_subleading_jet_eta_aa_ttbar, hist_leading_jet_eta_aa_ttbar, hist_jet_centrality_aa_ttbar,
 hist_delta_eta_jj_aa_ttbar, hist_m_w_leptonic_aa_ttbar, hist_m_w_hadronic_aa_ttbar,
 background_efficiency_pre_aa_ttbar, background_efficiency_final_aa_ttbar) = process_file(
    aa_ttbar_background_file_path, hist_lepton_aa_ttbar, hist_leading_jet_aa_ttbar, hist_lepton_eta_aa_ttbar,
    hist_delta_r_aa_ttbar, hist_missing_et_aa_ttbar, hist_subleading_jet_eta_aa_ttbar, hist_leading_jet_eta_aa_ttbar,
    hist_jet_centrality_aa_ttbar, hist_delta_eta_jj_aa_ttbar, hist_m_w_leptonic_aa_ttbar, hist_m_w_hadronic_aa_ttbar
)

# === Process the AA_TAUTAU background file ===
(hist_lepton_aa_tautau, hist_leading_jet_aa_tautau, hist_lepton_eta_aa_tautau, hist_delta_r_aa_tautau,
 hist_missing_et_aa_tautau, hist_subleading_jet_eta_aa_tautau, hist_leading_jet_eta_aa_tautau, hist_jet_centrality_aa_tautau,
 hist_delta_eta_jj_aa_tautau, hist_m_w_leptonic_aa_tautau, hist_m_w_hadronic_aa_tautau,
 background_efficiency_pre_aa_tautau, background_efficiency_final_aa_tautau) = process_file(
    aa_tautau_background_file_path, hist_lepton_aa_tautau, hist_leading_jet_aa_tautau, hist_lepton_eta_aa_tautau,
    hist_delta_r_aa_tautau, hist_missing_et_aa_tautau, hist_subleading_jet_eta_aa_tautau, hist_leading_jet_eta_aa_tautau,
    hist_jet_centrality_aa_tautau, hist_delta_eta_jj_aa_tautau, hist_m_w_leptonic_aa_tautau, hist_m_w_hadronic_aa_tautau
)

# === Process the AA_MUMU background file ===
(hist_lepton_aa_mumu, hist_leading_jet_aa_mumu, hist_lepton_eta_aa_mumu, hist_delta_r_aa_mumu,
 hist_missing_et_aa_mumu, hist_subleading_jet_eta_aa_mumu, hist_leading_jet_eta_aa_mumu, hist_jet_centrality_aa_mumu,
 hist_delta_eta_jj_aa_mumu, hist_m_w_leptonic_aa_mumu, hist_m_w_hadronic_aa_mumu,
 background_efficiency_pre_aa_mumu, background_efficiency_final_aa_mumu) = process_file(
    aa_mumu_background_file_path, hist_lepton_aa_mumu, hist_leading_jet_aa_mumu, hist_lepton_eta_aa_mumu,
    hist_delta_r_aa_mumu, hist_missing_et_aa_mumu, hist_subleading_jet_eta_aa_mumu, hist_leading_jet_eta_aa_mumu,
    hist_jet_centrality_aa_mumu, hist_delta_eta_jj_aa_mumu, hist_m_w_leptonic_aa_mumu, hist_m_w_hadronic_aa_mumu
)

# === Process the INCLUSIVE_TTBAR background file ===
(hist_lepton_inclusive_ttbar, hist_leading_jet_inclusive_ttbar, hist_lepton_eta_inclusive_ttbar, hist_delta_r_inclusive_ttbar,
 hist_missing_et_inclusive_ttbar, hist_subleading_jet_eta_inclusive_ttbar, hist_leading_jet_eta_inclusive_ttbar, hist_jet_centrality_inclusive_ttbar,
 hist_delta_eta_jj_inclusive_ttbar, hist_m_w_leptonic_inclusive_ttbar, hist_m_w_hadronic_inclusive_ttbar,
 background_efficiency_pre_inclusive_ttbar, background_efficiency_final_inclusive_ttbar) = process_file(
    inclusive_ttbar_background_file_path, hist_lepton_inclusive_ttbar, hist_leading_jet_inclusive_ttbar, hist_lepton_eta_inclusive_ttbar,
    hist_delta_r_inclusive_ttbar, hist_missing_et_inclusive_ttbar, hist_subleading_jet_eta_inclusive_ttbar, hist_leading_jet_eta_inclusive_ttbar,
    hist_jet_centrality_inclusive_ttbar, hist_delta_eta_jj_inclusive_ttbar, hist_m_w_leptonic_inclusive_ttbar, hist_m_w_hadronic_inclusive_ttbar
)

# === Process the SINGLE_TOP background file ===
(hist_lepton_single_top, hist_leading_jet_single_top, hist_lepton_eta_single_top, hist_delta_r_single_top,
 hist_missing_et_single_top, hist_subleading_jet_eta_single_top, hist_leading_jet_eta_single_top, hist_jet_centrality_single_top,
 hist_delta_eta_jj_single_top, hist_m_w_leptonic_single_top, hist_m_w_hadronic_single_top,
 background_efficiency_pre_single_top, background_efficiency_final_single_top) = process_file(
    single_top_background_file_path, hist_lepton_single_top, hist_leading_jet_single_top, hist_lepton_eta_single_top,
    hist_delta_r_single_top, hist_missing_et_single_top, hist_subleading_jet_eta_single_top, hist_leading_jet_eta_single_top,
    hist_jet_centrality_single_top, hist_delta_eta_jj_single_top, hist_m_w_leptonic_single_top, hist_m_w_hadronic_single_top
)

# === Process the W_PRODUCTION background file ===
(hist_lepton_w_production, hist_leading_jet_w_production, hist_lepton_eta_w_production, hist_delta_r_w_production,
 hist_missing_et_w_production, hist_subleading_jet_eta_w_production, hist_leading_jet_eta_w_production, hist_jet_centrality_w_production,
 hist_delta_eta_jj_w_production, hist_m_w_leptonic_w_production, hist_m_w_hadronic_w_production,
 background_efficiency_pre_w_production, background_efficiency_final_w_production) = process_file(
    w_production_background_file_path, hist_lepton_w_production, hist_leading_jet_w_production, hist_lepton_eta_w_production,
    hist_delta_r_w_production, hist_missing_et_w_production, hist_subleading_jet_eta_w_production, hist_leading_jet_eta_w_production,
    hist_jet_centrality_w_production, hist_delta_eta_jj_w_production, hist_m_w_leptonic_w_production, hist_m_w_hadronic_w_production
)

# === Process the Z_PRODUCTION background file ===
(hist_lepton_z_production, hist_leading_jet_z_production, hist_lepton_eta_z_production, hist_delta_r_z_production,
 hist_missing_et_z_production, hist_subleading_jet_eta_z_production, hist_leading_jet_eta_z_production, hist_jet_centrality_z_production,
 hist_delta_eta_jj_z_production, hist_m_w_leptonic_z_production, hist_m_w_hadronic_z_production,
 background_efficiency_pre_z_production, background_efficiency_final_z_production) = process_file(
    z_production_background_file_path, hist_lepton_z_production, hist_leading_jet_z_production, hist_lepton_eta_z_production,
    hist_delta_r_z_production, hist_missing_et_z_production, hist_subleading_jet_eta_z_production, hist_leading_jet_eta_z_production,
    hist_jet_centrality_z_production, hist_delta_eta_jj_z_production, hist_m_w_leptonic_z_production, hist_m_w_hadronic_z_production
)


# === Process the WWJ_PRODUCTION background file ===

(hist_lepton_wwj_production, hist_leading_jet_wwj_production, hist_lepton_eta_wwj_production, hist_delta_r_wwj_production,
 hist_missing_et_wwj_production, hist_subleading_jet_eta_wwj_production, hist_leading_jet_eta_wwj_production, hist_jet_centrality_wwj_production,
 hist_delta_eta_jj_wwj_production, hist_m_w_leptonic_wwj_production, hist_m_w_hadronic_wwj_production,
 background_efficiency_pre_wwj_production, background_efficiency_final_wwj_production) = process_file(
    wwj_production_background_file_path, hist_lepton_wwj_production, hist_leading_jet_wwj_production, hist_lepton_eta_wwj_production,
    hist_delta_r_wwj_production, hist_missing_et_wwj_production, hist_subleading_jet_eta_wwj_production, hist_leading_jet_eta_wwj_production,
    hist_jet_centrality_wwj_production, hist_delta_eta_jj_wwj_production, hist_m_w_leptonic_wwj_production, hist_m_w_hadronic_wwj_production
)


# === Process the ZZJ_PRODUCTION background file ===

(hist_lepton_zzj_production, hist_leading_jet_zzj_production, hist_lepton_eta_zzj_production, hist_delta_r_zzj_production,
 hist_missing_et_zzj_production, hist_subleading_jet_eta_zzj_production, hist_leading_jet_eta_zzj_production, hist_jet_centrality_zzj_production,
 hist_delta_eta_jj_zzj_production, hist_m_w_leptonic_zzj_production, hist_m_w_hadronic_zzj_production,
 background_efficiency_pre_zzj_production, background_efficiency_final_zzj_production) = process_file(
    zzj_production_background_file_path, hist_lepton_zzj_production, hist_leading_jet_zzj_production, hist_lepton_eta_zzj_production,
    hist_delta_r_zzj_production, hist_missing_et_zzj_production, hist_subleading_jet_eta_zzj_production, hist_leading_jet_eta_zzj_production,
    hist_jet_centrality_zzj_production, hist_delta_eta_jj_zzj_production, hist_m_w_leptonic_zzj_production, hist_m_w_hadronic_zzj_production
)



# === Process the WZJ_PRODUCTION background file ===

(hist_lepton_wzj_production, hist_leading_jet_wzj_production, hist_lepton_eta_wzj_production, hist_delta_r_wzj_production,
 hist_missing_et_wzj_production, hist_subleading_jet_eta_wzj_production, hist_leading_jet_eta_wzj_production, hist_jet_centrality_wzj_production,
 hist_delta_eta_jj_wzj_production, hist_m_w_leptonic_wzj_production, hist_m_w_hadronic_wzj_production,
 background_efficiency_pre_wzj_production, background_efficiency_final_wzj_production) = process_file(
    wzj_production_background_file_path, hist_lepton_wzj_production, hist_leading_jet_wzj_production, hist_lepton_eta_wzj_production,
    hist_delta_r_wzj_production, hist_missing_et_wzj_production, hist_subleading_jet_eta_wzj_production, hist_leading_jet_eta_wzj_production,
    hist_jet_centrality_wzj_production, hist_delta_eta_jj_wzj_production, hist_m_w_leptonic_wzj_production, hist_m_w_hadronic_wzj_production
)









# Print selection efficiencies for all signals
print("\n=== Signal Selection Efficiencies ===")
for signal_name, efficiencies in signal_efficiencies.items():
    print(f"{signal_name} Selection Efficiency (Pre)  : {efficiencies['efficiency_pre']:.2%}")
    print(f"{signal_name} Selection Efficiency (Final): {efficiencies['efficiency_final']:.2%}")
    print("-" * 50)





# Print background efficiencies
print("\n=== Background Selection Efficiencies ===")
print(f"aa_ww          : Pre = {background_efficiency_pre_aa_ww:.2%}, Final = {background_efficiency_final_aa_ww:.2%}")
print(f"aa_ttbar       : Pre = {background_efficiency_pre_aa_ttbar:.2%}, Final = {background_efficiency_final_aa_ttbar:.2%}")
print(f"aa_tautau      : Pre = {background_efficiency_pre_aa_tautau:.2%}, Final = {background_efficiency_final_aa_tautau:.2%}")
print(f"aa_mumu        : Pre = {background_efficiency_pre_aa_mumu:.2%}, Final = {background_efficiency_final_aa_mumu:.2%}")
print(f"inclusive_ttbar: Pre = {background_efficiency_pre_inclusive_ttbar:.2%}, Final = {background_efficiency_final_inclusive_ttbar:.2%}")
print(f"single_top     : Pre = {background_efficiency_pre_single_top:.2%}, Final = {background_efficiency_final_single_top:.2%}")
print(f"w_production   : Pre = {background_efficiency_pre_w_production:.2%}, Final = {background_efficiency_final_w_production:.2%}")
print(f"z_production   : Pre = {background_efficiency_pre_z_production:.2%}, Final = {background_efficiency_final_z_production:.2%}")
print(f"wwj_production : Pre = {background_efficiency_pre_wwj_production:.2%}, Final = {background_efficiency_final_wwj_production:.2%}")
print(f"zzj_production : Pre = {background_efficiency_pre_zzj_production:.2%}, Final = {background_efficiency_final_zzj_production:.2%}")
print(f"wzj_production : Pre = {background_efficiency_pre_wzj_production:.2%}, Final = {background_efficiency_final_wzj_production:.2%}")
print("=" * 50)






# Dictionary to store differential cross-sections for each signal
signal_dsigma = {}


# Calculate differential cross-sections for all signals
for signal_name, histograms in signal_histograms.items():
    print(f"Calculating differential cross-sections for {signal_name}...")

    cross_section = signal_cross_sections[signal_name]  # Get cross-section for this signal

    signal_dsigma[signal_name] = {
        "pt_bins_lepton": calculate_dsigma(histograms["hist_lepton_pt"], cross_section, bin_width_pt_lepton),
        "pt_bins_jet": calculate_dsigma(histograms["hist_leading_jet_pt"], cross_section, bin_width_pt_jet),
        "eta_bins_lepton": calculate_dsigma(histograms["hist_lepton_eta"], cross_section, bin_width_eta),
        "delta_r_bins": calculate_dsigma(histograms["hist_delta_r"], cross_section, bin_width_delta_r),
        "met_bins": calculate_dsigma(histograms["hist_missing_et"], cross_section, bin_width_met),
        "centrality_bins": calculate_dsigma(histograms["hist_subleading_jet_eta"], cross_section, bin_width_centrality),
        "exp_centrality_bins": calculate_dsigma(histograms["hist_leading_jet_eta"], cross_section, bin_width_exp_centrality),
        "jet_centrality_bins": calculate_dsigma(histograms["hist_jet_centrality"], cross_section, bin_width_jet_centrality),
        "delta_eta_jj_bins": calculate_dsigma(histograms["hist_delta_eta_jj"], cross_section, bin_width_delta_eta_jj),
        "m_w_hadronic_bins": calculate_dsigma(histograms["hist_m_w_hadronic"], cross_section, bin_width_m_w_hadronic),
        "m_w_leptonic_bins": calculate_dsigma(histograms["hist_m_w_leptonic"], cross_section, bin_width_m_w_leptonic),
    }





# === Differential cross-sections for aa_ww background ===
background_dsigma_aa_ww = {
    "pt_bins_lepton": calculate_dsigma(hist_lepton_aa_ww, aa_ww_background_cross_section, bin_width_pt_lepton),
    "pt_bins_jet": calculate_dsigma(hist_leading_jet_aa_ww, aa_ww_background_cross_section, bin_width_pt_jet),
    "eta_bins_lepton": calculate_dsigma(hist_lepton_eta_aa_ww, aa_ww_background_cross_section, bin_width_eta),
    "delta_r_bins": calculate_dsigma(hist_delta_r_aa_ww, aa_ww_background_cross_section, bin_width_delta_r),
    "met_bins": calculate_dsigma(hist_missing_et_aa_ww, aa_ww_background_cross_section, bin_width_met),
    "centrality_bins": calculate_dsigma(hist_subleading_jet_eta_aa_ww, aa_ww_background_cross_section, bin_width_centrality),
    "exp_centrality_bins": calculate_dsigma(hist_leading_jet_eta_aa_ww, aa_ww_background_cross_section, bin_width_exp_centrality),
    "jet_centrality_bins": calculate_dsigma(hist_jet_centrality_aa_ww, aa_ww_background_cross_section, bin_width_jet_centrality),
    "delta_eta_jj_bins": calculate_dsigma(hist_delta_eta_jj_aa_ww, aa_ww_background_cross_section, bin_width_delta_eta_jj),
    "m_w_hadronic_bins": calculate_dsigma(hist_m_w_hadronic_aa_ww, aa_ww_background_cross_section, bin_width_m_w_hadronic),
    "m_w_leptonic_bins": calculate_dsigma(hist_m_w_leptonic_aa_ww, aa_ww_background_cross_section, bin_width_m_w_leptonic),
}

# === Differential cross-sections for aa_ttbar background ===
background_dsigma_aa_ttbar = {
    "pt_bins_lepton": calculate_dsigma(hist_lepton_aa_ttbar, aa_ttbar_background_cross_section, bin_width_pt_lepton),
    "pt_bins_jet": calculate_dsigma(hist_leading_jet_aa_ttbar, aa_ttbar_background_cross_section, bin_width_pt_jet),
    "eta_bins_lepton": calculate_dsigma(hist_lepton_eta_aa_ttbar, aa_ttbar_background_cross_section, bin_width_eta),
    "delta_r_bins": calculate_dsigma(hist_delta_r_aa_ttbar, aa_ttbar_background_cross_section, bin_width_delta_r),
    "met_bins": calculate_dsigma(hist_missing_et_aa_ttbar, aa_ttbar_background_cross_section, bin_width_met),
    "centrality_bins": calculate_dsigma(hist_subleading_jet_eta_aa_ttbar, aa_ttbar_background_cross_section, bin_width_centrality),
    "exp_centrality_bins": calculate_dsigma(hist_leading_jet_eta_aa_ttbar, aa_ttbar_background_cross_section, bin_width_exp_centrality),
    "jet_centrality_bins": calculate_dsigma(hist_jet_centrality_aa_ttbar, aa_ttbar_background_cross_section, bin_width_jet_centrality),
    "delta_eta_jj_bins": calculate_dsigma(hist_delta_eta_jj_aa_ttbar, aa_ttbar_background_cross_section, bin_width_delta_eta_jj),
    "m_w_hadronic_bins": calculate_dsigma(hist_m_w_hadronic_aa_ttbar, aa_ttbar_background_cross_section, bin_width_m_w_hadronic),
    "m_w_leptonic_bins": calculate_dsigma(hist_m_w_leptonic_aa_ttbar, aa_ttbar_background_cross_section, bin_width_m_w_leptonic),
}

# === Differential cross-sections for aa_tautau background ===
background_dsigma_aa_tautau = {
    "pt_bins_lepton": calculate_dsigma(hist_lepton_aa_tautau, aa_tautau_background_cross_section, bin_width_pt_lepton),
    "pt_bins_jet": calculate_dsigma(hist_leading_jet_aa_tautau, aa_tautau_background_cross_section, bin_width_pt_jet),
    "eta_bins_lepton": calculate_dsigma(hist_lepton_eta_aa_tautau, aa_tautau_background_cross_section, bin_width_eta),
    "delta_r_bins": calculate_dsigma(hist_delta_r_aa_tautau, aa_tautau_background_cross_section, bin_width_delta_r),
    "met_bins": calculate_dsigma(hist_missing_et_aa_tautau, aa_tautau_background_cross_section, bin_width_met),
    "centrality_bins": calculate_dsigma(hist_subleading_jet_eta_aa_tautau, aa_tautau_background_cross_section, bin_width_centrality),
    "exp_centrality_bins": calculate_dsigma(hist_leading_jet_eta_aa_tautau, aa_tautau_background_cross_section, bin_width_exp_centrality),
    "jet_centrality_bins": calculate_dsigma(hist_jet_centrality_aa_tautau, aa_tautau_background_cross_section, bin_width_jet_centrality),
    "delta_eta_jj_bins": calculate_dsigma(hist_delta_eta_jj_aa_tautau, aa_tautau_background_cross_section, bin_width_delta_eta_jj),
    "m_w_hadronic_bins": calculate_dsigma(hist_m_w_hadronic_aa_tautau, aa_tautau_background_cross_section, bin_width_m_w_hadronic),
    "m_w_leptonic_bins": calculate_dsigma(hist_m_w_leptonic_aa_tautau, aa_tautau_background_cross_section, bin_width_m_w_leptonic),
}

# === Differential cross-sections for aa_mumu background ===
background_dsigma_aa_mumu = {
    "pt_bins_lepton": calculate_dsigma(hist_lepton_aa_mumu, aa_mumu_background_cross_section, bin_width_pt_lepton),
    "pt_bins_jet": calculate_dsigma(hist_leading_jet_aa_mumu, aa_mumu_background_cross_section, bin_width_pt_jet),
    "eta_bins_lepton": calculate_dsigma(hist_lepton_eta_aa_mumu, aa_mumu_background_cross_section, bin_width_eta),
    "delta_r_bins": calculate_dsigma(hist_delta_r_aa_mumu, aa_mumu_background_cross_section, bin_width_delta_r),
    "met_bins": calculate_dsigma(hist_missing_et_aa_mumu, aa_mumu_background_cross_section, bin_width_met),
    "centrality_bins": calculate_dsigma(hist_subleading_jet_eta_aa_mumu, aa_mumu_background_cross_section, bin_width_centrality),
    "exp_centrality_bins": calculate_dsigma(hist_leading_jet_eta_aa_mumu, aa_mumu_background_cross_section, bin_width_exp_centrality),
    "jet_centrality_bins": calculate_dsigma(hist_jet_centrality_aa_mumu, aa_mumu_background_cross_section, bin_width_jet_centrality),
    "delta_eta_jj_bins": calculate_dsigma(hist_delta_eta_jj_aa_mumu, aa_mumu_background_cross_section, bin_width_delta_eta_jj),
    "m_w_hadronic_bins": calculate_dsigma(hist_m_w_hadronic_aa_mumu, aa_mumu_background_cross_section, bin_width_m_w_hadronic),
    "m_w_leptonic_bins": calculate_dsigma(hist_m_w_leptonic_aa_mumu, aa_mumu_background_cross_section, bin_width_m_w_leptonic),
}

# === Differential cross-sections for inclusive_ttbar background ===
background_dsigma_inclusive_ttbar = {
    "pt_bins_lepton": calculate_dsigma(hist_lepton_inclusive_ttbar, inclusive_ttbar_background_cross_section, bin_width_pt_lepton),
    "pt_bins_jet": calculate_dsigma(hist_leading_jet_inclusive_ttbar, inclusive_ttbar_background_cross_section, bin_width_pt_jet),
    "eta_bins_lepton": calculate_dsigma(hist_lepton_eta_inclusive_ttbar, inclusive_ttbar_background_cross_section, bin_width_eta),
    "delta_r_bins": calculate_dsigma(hist_delta_r_inclusive_ttbar, inclusive_ttbar_background_cross_section, bin_width_delta_r),
    "met_bins": calculate_dsigma(hist_missing_et_inclusive_ttbar, inclusive_ttbar_background_cross_section, bin_width_met),
    "centrality_bins": calculate_dsigma(hist_subleading_jet_eta_inclusive_ttbar, inclusive_ttbar_background_cross_section, bin_width_centrality),
    "exp_centrality_bins": calculate_dsigma(hist_leading_jet_eta_inclusive_ttbar, inclusive_ttbar_background_cross_section, bin_width_exp_centrality),
    "jet_centrality_bins": calculate_dsigma(hist_jet_centrality_inclusive_ttbar, inclusive_ttbar_background_cross_section, bin_width_jet_centrality),
    "delta_eta_jj_bins": calculate_dsigma(hist_delta_eta_jj_inclusive_ttbar, inclusive_ttbar_background_cross_section, bin_width_delta_eta_jj),
    "m_w_hadronic_bins": calculate_dsigma(hist_m_w_hadronic_inclusive_ttbar, inclusive_ttbar_background_cross_section, bin_width_m_w_hadronic),
    "m_w_leptonic_bins": calculate_dsigma(hist_m_w_leptonic_inclusive_ttbar, inclusive_ttbar_background_cross_section, bin_width_m_w_leptonic),
}

# === Differential cross-sections for single_top background ===
background_dsigma_single_top = {
    "pt_bins_lepton": calculate_dsigma(hist_lepton_single_top, single_top_background_cross_section, bin_width_pt_lepton),
    "pt_bins_jet": calculate_dsigma(hist_leading_jet_single_top, single_top_background_cross_section, bin_width_pt_jet),
    "eta_bins_lepton": calculate_dsigma(hist_lepton_eta_single_top, single_top_background_cross_section, bin_width_eta),
    "delta_r_bins": calculate_dsigma(hist_delta_r_single_top, single_top_background_cross_section, bin_width_delta_r),
    "met_bins": calculate_dsigma(hist_missing_et_single_top, single_top_background_cross_section, bin_width_met),
    "centrality_bins": calculate_dsigma(hist_subleading_jet_eta_single_top, single_top_background_cross_section, bin_width_centrality),
    "exp_centrality_bins": calculate_dsigma(hist_leading_jet_eta_single_top, single_top_background_cross_section, bin_width_exp_centrality),
    "jet_centrality_bins": calculate_dsigma(hist_jet_centrality_single_top, single_top_background_cross_section, bin_width_jet_centrality),
    "delta_eta_jj_bins": calculate_dsigma(hist_delta_eta_jj_single_top, single_top_background_cross_section, bin_width_delta_eta_jj),
    "m_w_hadronic_bins": calculate_dsigma(hist_m_w_hadronic_single_top, single_top_background_cross_section, bin_width_m_w_hadronic),
    "m_w_leptonic_bins": calculate_dsigma(hist_m_w_leptonic_single_top, single_top_background_cross_section, bin_width_m_w_leptonic),
}

# === Differential cross-sections for w_production background ===
background_dsigma_w_production = {
    "pt_bins_lepton": calculate_dsigma(hist_lepton_w_production, w_production_background_cross_section, bin_width_pt_lepton),
    "pt_bins_jet": calculate_dsigma(hist_leading_jet_w_production, w_production_background_cross_section, bin_width_pt_jet),
    "eta_bins_lepton": calculate_dsigma(hist_lepton_eta_w_production, w_production_background_cross_section, bin_width_eta),
    "delta_r_bins": calculate_dsigma(hist_delta_r_w_production, w_production_background_cross_section, bin_width_delta_r),
    "met_bins": calculate_dsigma(hist_missing_et_w_production, w_production_background_cross_section, bin_width_met),
    "centrality_bins": calculate_dsigma(hist_subleading_jet_eta_w_production, w_production_background_cross_section, bin_width_centrality),
    "exp_centrality_bins": calculate_dsigma(hist_leading_jet_eta_w_production, w_production_background_cross_section, bin_width_exp_centrality),
    "jet_centrality_bins": calculate_dsigma(hist_jet_centrality_w_production, w_production_background_cross_section, bin_width_jet_centrality),
    "delta_eta_jj_bins": calculate_dsigma(hist_delta_eta_jj_w_production, w_production_background_cross_section, bin_width_delta_eta_jj),
    "m_w_hadronic_bins": calculate_dsigma(hist_m_w_hadronic_w_production, w_production_background_cross_section, bin_width_m_w_hadronic),
    "m_w_leptonic_bins": calculate_dsigma(hist_m_w_leptonic_w_production, w_production_background_cross_section, bin_width_m_w_leptonic),
}

# === Differential cross-sections for z_production background ===
background_dsigma_z_production = {
    "pt_bins_lepton": calculate_dsigma(hist_lepton_z_production, z_production_background_cross_section, bin_width_pt_lepton),
    "pt_bins_jet": calculate_dsigma(hist_leading_jet_z_production, z_production_background_cross_section, bin_width_pt_jet),
    "eta_bins_lepton": calculate_dsigma(hist_lepton_eta_z_production, z_production_background_cross_section, bin_width_eta),
    "delta_r_bins": calculate_dsigma(hist_delta_r_z_production, z_production_background_cross_section, bin_width_delta_r),
    "met_bins": calculate_dsigma(hist_missing_et_z_production, z_production_background_cross_section, bin_width_met),
    "centrality_bins": calculate_dsigma(hist_subleading_jet_eta_z_production, z_production_background_cross_section, bin_width_centrality),
    "exp_centrality_bins": calculate_dsigma(hist_leading_jet_eta_z_production, z_production_background_cross_section, bin_width_exp_centrality),
    "jet_centrality_bins": calculate_dsigma(hist_jet_centrality_z_production, z_production_background_cross_section, bin_width_jet_centrality),
    "delta_eta_jj_bins": calculate_dsigma(hist_delta_eta_jj_z_production, z_production_background_cross_section, bin_width_delta_eta_jj),
    "m_w_hadronic_bins": calculate_dsigma(hist_m_w_hadronic_z_production, z_production_background_cross_section, bin_width_m_w_hadronic),
    "m_w_leptonic_bins": calculate_dsigma(hist_m_w_leptonic_z_production, z_production_background_cross_section, bin_width_m_w_leptonic),
}




# === Differential cross-sections for wwj_production background ===
background_dsigma_wwj_production = {
    "pt_bins_lepton": calculate_dsigma(hist_lepton_wwj_production, wwj_production_background_cross_section, bin_width_pt_lepton),
    "pt_bins_jet": calculate_dsigma(hist_leading_jet_wwj_production, wwj_production_background_cross_section, bin_width_pt_jet),
    "eta_bins_lepton": calculate_dsigma(hist_lepton_eta_wwj_production, wwj_production_background_cross_section, bin_width_eta),
    "delta_r_bins": calculate_dsigma(hist_delta_r_wwj_production, wwj_production_background_cross_section, bin_width_delta_r),
    "met_bins": calculate_dsigma(hist_missing_et_wwj_production, wwj_production_background_cross_section, bin_width_met),
    "centrality_bins": calculate_dsigma(hist_subleading_jet_eta_wwj_production, wwj_production_background_cross_section, bin_width_centrality),
    "exp_centrality_bins": calculate_dsigma(hist_leading_jet_eta_wwj_production, wwj_production_background_cross_section, bin_width_exp_centrality),
    "jet_centrality_bins": calculate_dsigma(hist_jet_centrality_wwj_production, wwj_production_background_cross_section, bin_width_jet_centrality),
    "delta_eta_jj_bins": calculate_dsigma(hist_delta_eta_jj_wwj_production, wwj_production_background_cross_section, bin_width_delta_eta_jj),
    "m_w_hadronic_bins": calculate_dsigma(hist_m_w_hadronic_wwj_production, wwj_production_background_cross_section, bin_width_m_w_hadronic),
    "m_w_leptonic_bins": calculate_dsigma(hist_m_w_leptonic_wwj_production, wwj_production_background_cross_section, bin_width_m_w_leptonic),
}




# === Differential cross-sections for zzj_production background ===
background_dsigma_zzj_production = {
    "pt_bins_lepton": calculate_dsigma(hist_lepton_zzj_production, zzj_production_background_cross_section, bin_width_pt_lepton),
    "pt_bins_jet": calculate_dsigma(hist_leading_jet_zzj_production, zzj_production_background_cross_section, bin_width_pt_jet),
    "eta_bins_lepton": calculate_dsigma(hist_lepton_eta_zzj_production, zzj_production_background_cross_section, bin_width_eta),
    "delta_r_bins": calculate_dsigma(hist_delta_r_zzj_production, zzj_production_background_cross_section, bin_width_delta_r),
    "met_bins": calculate_dsigma(hist_missing_et_zzj_production, zzj_production_background_cross_section, bin_width_met),
    "centrality_bins": calculate_dsigma(hist_subleading_jet_eta_zzj_production, zzj_production_background_cross_section, bin_width_centrality),
    "exp_centrality_bins": calculate_dsigma(hist_leading_jet_eta_zzj_production, zzj_production_background_cross_section, bin_width_exp_centrality),
    "jet_centrality_bins": calculate_dsigma(hist_jet_centrality_zzj_production, zzj_production_background_cross_section, bin_width_jet_centrality),
    "delta_eta_jj_bins": calculate_dsigma(hist_delta_eta_jj_zzj_production, zzj_production_background_cross_section, bin_width_delta_eta_jj),
    "m_w_hadronic_bins": calculate_dsigma(hist_m_w_hadronic_zzj_production, zzj_production_background_cross_section, bin_width_m_w_hadronic),
    "m_w_leptonic_bins": calculate_dsigma(hist_m_w_leptonic_zzj_production, zzj_production_background_cross_section, bin_width_m_w_leptonic),
}



# === Differential cross-sections for wzj_production background ===
background_dsigma_wzj_production = {
    "pt_bins_lepton": calculate_dsigma(hist_lepton_wzj_production, wzj_production_background_cross_section, bin_width_pt_lepton),
    "pt_bins_jet": calculate_dsigma(hist_leading_jet_wzj_production, wzj_production_background_cross_section, bin_width_pt_jet),
    "eta_bins_lepton": calculate_dsigma(hist_lepton_eta_wzj_production, wzj_production_background_cross_section, bin_width_eta),
    "delta_r_bins": calculate_dsigma(hist_delta_r_wzj_production, wzj_production_background_cross_section, bin_width_delta_r),
    "met_bins": calculate_dsigma(hist_missing_et_wzj_production, wzj_production_background_cross_section, bin_width_met),
    "centrality_bins": calculate_dsigma(hist_subleading_jet_eta_wzj_production, wzj_production_background_cross_section, bin_width_centrality),
    "exp_centrality_bins": calculate_dsigma(hist_leading_jet_eta_wzj_production, wzj_production_background_cross_section, bin_width_exp_centrality),
    "jet_centrality_bins": calculate_dsigma(hist_jet_centrality_wzj_production, wzj_production_background_cross_section, bin_width_jet_centrality),
    "delta_eta_jj_bins": calculate_dsigma(hist_delta_eta_jj_wzj_production, wzj_production_background_cross_section, bin_width_delta_eta_jj),
    "m_w_hadronic_bins": calculate_dsigma(hist_m_w_hadronic_wzj_production, wzj_production_background_cross_section, bin_width_m_w_hadronic),
    "m_w_leptonic_bins": calculate_dsigma(hist_m_w_leptonic_wzj_production, wzj_production_background_cross_section, bin_width_m_w_leptonic),
}







#=========================================================================
#=========================================================================




plt.figure(figsize=(11, 12))  # Create a new figure for the leading jet η plot
plt.subplots_adjust(left=0.15, right=0.95, bottom=0.12, top=0.95)



# Define colors for each signal
signal_colors = {
    "$FM_{0} / \Lambda^4$": "green",
    "$FM_{1} / \Lambda^4$": "purple",
    "$FM_{2} / \Lambda^4$": "red",
    "$FM_{3} / \Lambda^4$": "orange"
}





# ===================================================
# ✅ Lepton \( p_T \) Differential Cross-Section (Semi-Leptonic)
# ===================================================



# ✅ Plot Signals
for signal_name, dsigma_data in signal_dsigma.items():
    pt_bins, dsigma = dsigma_data["pt_bins_lepton"]
    plt.step(pt_bins, dsigma, where="mid", alpha=0.7,
             label=f"Signal ($W^+W^-$) [{signal_name}]",
             color=signal_colors.get(signal_name, "black"), linewidth=3)

# ✅ Plot Backgrounds
# Each entry: (label, dsigma_dict, color, linestyle)
background_styles = [
    (r"SM ($\gamma\gamma \to W^+W^-$)", background_dsigma_aa_ww, "blue", "-"),
    (r"$\gamma\gamma \to t\bar{t}$ $(\times 10^2)$", background_dsigma_aa_ttbar, "teal", "--"),
    (r"$\gamma\gamma \to \tau^+\tau^-$", background_dsigma_aa_tautau, "orange", "-"),
    (r"$\gamma\gamma \to \mu^+\mu^-$", background_dsigma_aa_mumu, "brown", (0, (3, 1, 1, 1))),
    (r"Inclusive $t\bar{t}$", background_dsigma_inclusive_ttbar, "orchid", (0, (3, 1, 1, 1))),
    (r"Single Top", background_dsigma_single_top, "teal", "--"),
    (r"W Production", background_dsigma_w_production, "darkgreen", "-."),
    (r"Z Production", background_dsigma_z_production, "darkred", ":"),
    (r"$WWj$", background_dsigma_wwj_production, "royalblue", "--"),
    (r"$ZZj$ $(\times 10^2)$", background_dsigma_zzj_production, "mediumvioletred", "-."),
    (r"$WZj$", background_dsigma_wzj_production, "slategray", ":")
]


for label, dsigma_dict, color, style in background_styles:
    pt_bins, dsigma = dsigma_dict["pt_bins_lepton"]
    if sum(dsigma) > 0:
        plt.step(pt_bins, dsigma, where="mid", alpha=0.7,
                 label=label, color=color, linestyle=style, linewidth=3)

# ✅ Axis labels and title
plt.xlabel(r"$p_T^{\ell} \ \mathrm{[GeV]}$")
plt.ylabel(r"$\frac{d\sigma}{dp_T^{\ell}} \ \mathrm{[pb/GeV]}$")
plt.title(r"Delphes simulation : $e^- p \to e^- W^+ W^- p \to e^- j j \ell \nu_{\ell} p$ : LHeC@1.2 TeV", fontsize=20)
plt.yscale("log")
plt.ylim(1e-5, 1e-1)

# ✅ Legend, Grid, and Save
plt.legend()
plt.grid(True, linestyle="--", alpha=0.6)
plt.tight_layout()
plt.savefig("/home/hamzeh-khanpour/Documents/GitHub/LHeC_Fast_Simulation/Pythia8_Delphes_semi_leptonic_allsignal_bkgs/differential_cross_section_lepton_pt_allFMsignal_allbkgs_all_decay_mode.pdf", dpi=600)

plt.show()








# ===================================================
# ✅ Leading Jet \( p_T \) Differential Cross-Section (Semi-Leptonic)
# ===================================================

plt.figure(figsize=(11, 12))
plt.subplots_adjust(left=0.15, right=0.95, bottom=0.12, top=0.95)

# ✅ Plot Signals
for signal_name, dsigma_data in signal_dsigma.items():
    pt_bins, dsigma = dsigma_data["pt_bins_jet"]
    plt.step(pt_bins, dsigma, where="mid", alpha=0.7,
             label=f"Signal ($W^+W^-$) [{signal_name}]",
             color=signal_colors.get(signal_name, "black"), linewidth=3)

# ✅ Plot Backgrounds
background_styles = [
    (r"SM ($\gamma\gamma \to W^+W^-$)", background_dsigma_aa_ww, "blue", "-"),
    (r"$\gamma\gamma \to t\bar{t}$ $(\times 10^2)$", background_dsigma_aa_ttbar, "teal", "--"),
    (r"$\gamma\gamma \to \tau^+\tau^-$", background_dsigma_aa_tautau, "orange", "-"),
    (r"$\gamma\gamma \to \mu^+\mu^-$", background_dsigma_aa_mumu, "brown", (0, (3, 1, 1, 1))),
    (r"Inclusive $t\bar{t}$", background_dsigma_inclusive_ttbar, "orchid", (0, (3, 1, 1, 1))),
    (r"Single Top", background_dsigma_single_top, "teal", "--"),
    (r"W Production", background_dsigma_w_production, "darkgreen", "-."),
    (r"Z Production", background_dsigma_z_production, "darkred", ":"),
    (r"$WWj$", background_dsigma_wwj_production, "royalblue", (0, (5, 2))),
    (r"$ZZj$ $(\times 10^2)$", background_dsigma_zzj_production, "firebrick", (0, (1, 1))),
    (r"$WZj$", background_dsigma_wzj_production, "mediumslateblue", (0, (3, 2, 1, 2)))
]


for label, dsigma_dict, color, style in background_styles:
    pt_bins, dsigma = dsigma_dict["pt_bins_jet"]
    if sum(dsigma) > 0:
        plt.step(pt_bins, dsigma, where="mid", alpha=0.7,
                 label=label, color=color, linestyle=style, linewidth=3)

# ✅ Axis Labels and Title
plt.xlabel(r"$p_T^{\mathrm{leading~jet}} \ \mathrm{[GeV]}$")
plt.ylabel(r"$\frac{d\sigma}{dp_T^{\mathrm{leading~jet}}} \ \mathrm{[pb/GeV]}$")
plt.title(r"Delphes simulation : $e^- p \to e^- W^+ W^- p \to e^- j j \ell \nu_{\ell} p$ : LHeC@1.2 \ TeV", fontsize=20)
plt.yscale("log")
plt.ylim(1e-5, 1e-1)

# ✅ Legend and Layout
plt.legend()
plt.grid(True, linestyle="--", alpha=0.6)
plt.tight_layout()

# ✅ Save and Show
plt.savefig("/home/hamzeh-khanpour/Documents/GitHub/LHeC_Fast_Simulation/Pythia8_Delphes_semi_leptonic_allsignal_bkgs/differential_cross_section_jet_pt_allFMsignal_allbkgs_all_decay_mode.pdf", dpi=600)
plt.show()




# ===================================================
# ✅ Lepton \( \eta \) Differential Cross-Section (Semi-Leptonic)
# ===================================================

plt.figure(figsize=(11, 12))
plt.subplots_adjust(left=0.15, right=0.95, bottom=0.12, top=0.95)

# ✅ Plot Signals
for signal_name, dsigma_data in signal_dsigma.items():
    eta_bins, dsigma = dsigma_data["eta_bins_lepton"]
    plt.step(eta_bins, dsigma, where="mid", alpha=0.7,
             label=f"Signal ($W^+W^-$) [{signal_name}]",
             color=signal_colors.get(signal_name, "black"), linewidth=3)

# ✅ Plot Backgrounds
background_styles = [
    (r"SM ($\gamma\gamma \to W^+W^-$)", background_dsigma_aa_ww, "blue", "-"),
    (r"$\gamma\gamma \to t\bar{t}$ $(\times 10^2)$", background_dsigma_aa_ttbar, "teal", "--"),
    (r"$\gamma\gamma \to \tau^+\tau^-$", background_dsigma_aa_tautau, "orange", "-"),
    (r"$\gamma\gamma \to \mu^+\mu^-$", background_dsigma_aa_mumu, "brown", (0, (3, 1, 1, 1))),
    (r"Inclusive $t\bar{t}$", background_dsigma_inclusive_ttbar, "orchid", (0, (3, 1, 1, 1))),
    (r"Single Top", background_dsigma_single_top, "teal", "--"),
    (r"W Production", background_dsigma_w_production, "darkgreen", "-."),
    (r"Z Production", background_dsigma_z_production, "darkred", ":"),
    (r"$WWj$", background_dsigma_wwj_production, "royalblue", (0, (5, 2))),
    (r"$ZZj$ $(\times 10^2)$", background_dsigma_zzj_production, "firebrick", (0, (1, 1))),
    (r"$WZj$", background_dsigma_wzj_production, "mediumslateblue", (0, (3, 2, 1, 2)))
]


for label, dsigma_dict, color, style in background_styles:
    eta_bins, dsigma = dsigma_dict["eta_bins_lepton"]
    if sum(dsigma) > 0:
        plt.step(eta_bins, dsigma, where="mid", alpha=0.7,
                 label=label, color=color, linestyle=style, linewidth=3)

# ✅ Axis Labels and Title
plt.xlabel(r"$\eta^{\ell}$")
plt.ylabel(r"$\frac{d\sigma}{d\eta^{\ell}} \ \mathrm{[pb]}$")
plt.title(r"Delphes simulation : $e^- p \to e^- W^+ W^- p \to e^- j j \ell \nu_{\ell} p$ : LHeC@1.2 TeV", fontsize=20)
plt.yscale("log")
plt.ylim(1e-4, 10e1)

# ✅ Legend and Layout
plt.legend()
plt.grid(True, linestyle="--", alpha=0.6)
plt.tight_layout()

# ✅ Save and Show
plt.savefig("/home/hamzeh-khanpour/Documents/GitHub/LHeC_Fast_Simulation/Pythia8_Delphes_semi_leptonic_allsignal_bkgs/differential_cross_section_lepton_eta_allFMsignal_allbkgs_all_decay_mode.pdf", dpi=600)
plt.show()






# ===================================================
# ✅ ΔR(ℓ, jet) Differential Cross-Section (Semi-Leptonic)
# ===================================================

plt.figure(figsize=(11, 12))
plt.subplots_adjust(left=0.15, right=0.95, bottom=0.12, top=0.95)

# ✅ Plot Signals
for signal_name, dsigma_data in signal_dsigma.items():
    delta_r_bins, dsigma = dsigma_data["delta_r_bins"]
    plt.step(delta_r_bins, dsigma, where="mid", alpha=0.7,
             label=f"Signal ($W^+W^-$) [{signal_name}]",
             color=signal_colors.get(signal_name, "black"), linewidth=3)

# ✅ Plot Backgrounds
background_styles = [
    (r"SM ($\gamma\gamma \to W^+W^-$)", background_dsigma_aa_ww, "blue", "-"),
    (r"$\gamma\gamma \to t\bar{t}$ $(\times 10^2)$", background_dsigma_aa_ttbar, "teal", "--"),
    (r"$\gamma\gamma \to \tau^+\tau^-$", background_dsigma_aa_tautau, "orange", "-"),
    (r"$\gamma\gamma \to \mu^+\mu^-$", background_dsigma_aa_mumu, "brown", (0, (3, 1, 1, 1))),
    (r"Inclusive $t\bar{t}$", background_dsigma_inclusive_ttbar, "orchid", (0, (3, 1, 1, 1))),
    (r"Single Top", background_dsigma_single_top, "teal", "--"),
    (r"W Production", background_dsigma_w_production, "darkgreen", "-."),
    (r"Z Production", background_dsigma_z_production, "darkred", ":"),
    (r"$WWj$", background_dsigma_wwj_production, "royalblue", (0, (5, 2))),
    (r"$ZZj$ $(\times 10^2)$", background_dsigma_zzj_production, "crimson", (0, (1, 1))),
    (r"$WZj$", background_dsigma_wzj_production, "slateblue", (0, (3, 2, 1, 2)))
]


for label, dsigma_dict, color, style in background_styles:
    delta_r_bins, dsigma = dsigma_dict["delta_r_bins"]
    if sum(dsigma) > 0:
        plt.step(delta_r_bins, dsigma, where="mid", alpha=0.7,
                 label=label, color=color, linestyle=style, linewidth=3)

# ✅ Axis Labels and Title
plt.xlabel(r"$\Delta R(\ell, \mathrm{leading~jet})$")
plt.ylabel(r"$\frac{d\sigma}{d\Delta R} \ \mathrm{[pb]}$")
plt.title(r"Delphes simulation : $e^- p \to e^- W^+ W^- p \to e^- j j \ell \nu_{\ell} p$ : LHeC@1.2 TeV", fontsize=20)
plt.yscale("log")
plt.ylim(1e-4, 1.0)

# ✅ Legend and Layout
plt.legend()
plt.grid(True, linestyle="--", alpha=0.6)
plt.tight_layout()

# ✅ Save and Show
plt.savefig("/home/hamzeh-khanpour/Documents/GitHub/LHeC_Fast_Simulation/Pythia8_Delphes_semi_leptonic_allsignal_bkgs/differential_cross_section_delta_r_allFMsignal_allbkgs_all_decay_mode.pdf", dpi=600)
plt.show()








# ===================================================
# ✅ MET Differential Cross-Section (Semi-Leptonic)
# ===================================================

plt.figure(figsize=(11, 12))
plt.subplots_adjust(left=0.15, right=0.95, bottom=0.12, top=0.95)

# ✅ Plot Signals
for signal_name, dsigma_data in signal_dsigma.items():
    met_bins, dsigma = dsigma_data["met_bins"]
    plt.step(met_bins, dsigma, where="mid", alpha=0.7,
             label=f"Signal ($W^+W^-$) [{signal_name}]",
             color=signal_colors.get(signal_name, "black"), linewidth=3)

# ✅ Plot Backgrounds
background_styles = [
    (r"SM ($\gamma\gamma \to W^+W^-$)", background_dsigma_aa_ww, "blue", "-"),
    (r"$\gamma\gamma \to t\bar{t}$ $(\times 10^2)$", background_dsigma_aa_ttbar, "teal", "--"),
    (r"$\gamma\gamma \to \tau^+\tau^-$", background_dsigma_aa_tautau, "orange", "-"),
    (r"$\gamma\gamma \to \mu^+\mu^-$", background_dsigma_aa_mumu, "brown", (0, (3, 1, 1, 1))),
    (r"Inclusive $t\bar{t}$", background_dsigma_inclusive_ttbar, "orchid", (0, (3, 1, 1, 1))),
    (r"Single Top", background_dsigma_single_top, "teal", "--"),
    (r"W Production", background_dsigma_w_production, "darkgreen", "-."),
    (r"Z Production", background_dsigma_z_production, "darkred", ":"),
    (r"$WWj$", background_dsigma_wwj_production, "royalblue", (0, (5, 2))),
    (r"$ZZj$ $(\times 10^2)$", background_dsigma_zzj_production, "crimson", (0, (1, 1))),
    (r"$WZj$", background_dsigma_wzj_production, "slateblue", (0, (3, 2, 1, 2)))
]


for label, dsigma_dict, color, style in background_styles:
    met_bins, dsigma = dsigma_dict["met_bins"]
    if sum(dsigma) > 0:
        plt.step(met_bins, dsigma, where="mid", alpha=0.7,
                 label=label, color=color, linestyle=style, linewidth=3)

# ✅ Axis Labels and Title
plt.xlabel(r"$\mathrm{MET} \ \mathrm{[GeV]}$")
plt.ylabel(r"$\frac{d\sigma}{d\mathrm{MET}} \ \mathrm{[pb/GeV]}$")
plt.title(r"Delphes simulation : $e^- p \to e^- W^+ W^- p \to e^- j j \ell \nu_{\ell} p$ : LHeC@1.2 TeV", fontsize=20)
plt.yscale("log")
plt.ylim(1e-5, 1e-1)

# ✅ Legend and Layout
plt.legend()
plt.grid(True, linestyle="--", alpha=0.6)
plt.tight_layout()

# ✅ Save and Show
plt.savefig("/home/hamzeh-khanpour/Documents/GitHub/LHeC_Fast_Simulation/Pythia8_Delphes_semi_leptonic_allsignal_bkgs/differential_cross_section_met_allFMsignal_allbkgs_all_decay_mode.pdf", dpi=600)
plt.show()








# ===================================================
# ✅ Normalized Hadronic \( M_W^{jj} \) Distribution (Semi-Leptonic)
# ===================================================

plt.figure(figsize=(11, 12))
plt.subplots_adjust(left=0.15, right=0.95, bottom=0.12, top=0.95)

# ✅ Plot Normalized Signal Distributions
for signal_name, dsigma_data in signal_dsigma.items():
    m_w_hadronic_bins, dsigma = dsigma_data["m_w_hadronic_bins"]
    if np.sum(dsigma) > 0:
        normalized_dsigma = dsigma / np.sum(dsigma)
        plt.step(m_w_hadronic_bins, normalized_dsigma, where="mid", alpha=0.7,
                 label=f"Signal ($W^+W^-$) [{signal_name}]",
                 color=signal_colors.get(signal_name, "black"), linewidth=3)

# ✅ Plot Normalized Backgrounds
background_styles = [
    (r"SM ($\gamma\gamma \to W^+W^-$)", background_dsigma_aa_ww, "blue", "-"),
    (r"$\gamma\gamma \to t\bar{t}$ $(\times 10^2)$", background_dsigma_aa_ttbar, "purple", "-"),
    (r"$\gamma\gamma \to \tau^+\tau^-$", background_dsigma_aa_tautau, "orange", "-"),
    (r"$\gamma\gamma \to \mu^+\mu^-$", background_dsigma_aa_mumu, "brown", (0, (3, 1, 1, 1))),
    (r"Inclusive $t\bar{t}$", background_dsigma_inclusive_ttbar, "orchid", (0, (3, 1, 1, 1))),
    (r"Single Top", background_dsigma_single_top, "teal", "--"),
    (r"W Production", background_dsigma_w_production, "darkgreen", "-."),
    (r"Z Production", background_dsigma_z_production, "darkred", ":"),
    (r"$WWj$", background_dsigma_wwj_production, "royalblue", (0, (5, 2))),
    (r"$ZZj$ $(\times 10^2)$", background_dsigma_zzj_production, "crimson", (0, (1, 1))),
    (r"$WZj$", background_dsigma_wzj_production, "slateblue", (0, (3, 2, 1, 2)))
]



for label, dsigma_dict, color, style in background_styles:
    m_w_bins, dsigma = dsigma_dict["m_w_hadronic_bins"]
    if np.sum(dsigma) > 0:
        normalized_dsigma = dsigma / np.sum(dsigma)
        plt.step(m_w_bins, normalized_dsigma, where="mid", alpha=0.7,
                 label=label, color=color, linestyle=style, linewidth=3)

# ✅ Axis Labels and Title
plt.xlabel(r"$M_W^{\mathrm{j_1j_2}} \ \mathrm{[GeV]}$")
plt.ylabel("Normalized Entries")
plt.title(r"Delphes simulation : $e^- p \to e^- W^+ W^- p \to e^- j j \ell \nu_{\ell} p$ : LHeC@1.2 TeV", fontsize=20)
plt.ylim(0.0, 0.4)

# ✅ Legend and Layout
plt.legend()
plt.grid(True, linestyle="--", alpha=0.6)
plt.tight_layout()

# ✅ Save and Show
plt.savefig("/home/hamzeh-khanpour/Documents/GitHub/LHeC_Fast_Simulation/Pythia8_Delphes_semi_leptonic_allsignal_bkgs/normalized_m_w_hadronic_allFMsignal_allbkgs_all_decay_mode.pdf", dpi=600)
plt.show()







# ===================================================
# ✅ Normalized Leptonic \( M_W^{\ell \nu} \) Distribution (Semi-Leptonic)
# ===================================================

plt.figure(figsize=(11, 12))
plt.subplots_adjust(left=0.15, right=0.95, bottom=0.12, top=0.95)

# ✅ Plot Normalized Signal Distributions
for signal_name, dsigma_data in signal_dsigma.items():
    m_w_leptonic_bins, dsigma = dsigma_data["m_w_leptonic_bins"]
    if np.sum(dsigma) > 0:
        normalized_dsigma = dsigma / np.sum(dsigma)
        plt.step(m_w_leptonic_bins, normalized_dsigma, where="mid", alpha=0.7,
                 label=f"Signal ($W^+W^-$) [{signal_name}]",
                 color=signal_colors.get(signal_name, "black"), linewidth=3)

# ✅ Plot Normalized Backgrounds
background_styles = [
    (r"SM ($\gamma\gamma \to W^+W^-$)", background_dsigma_aa_ww, "blue", "-"),
    (r"$\gamma\gamma \to t\bar{t}$ $(\times 10^2)$", background_dsigma_aa_ttbar, "teal", "--"),
    (r"$\gamma\gamma \to \tau^+\tau^-$", background_dsigma_aa_tautau, "orange", "-"),
    (r"$\gamma\gamma \to \mu^+\mu^-$", background_dsigma_aa_mumu, "brown", (0, (3, 1, 1, 1))),
    (r"Inclusive $t\bar{t}$", background_dsigma_inclusive_ttbar, "orchid", (0, (3, 1, 1, 1))),
    (r"Single Top", background_dsigma_single_top, "teal", "--"),
    (r"W Production", background_dsigma_w_production, "darkgreen", "-."),
    (r"Z Production", background_dsigma_z_production, "darkred", ":"),
    (r"$WWj$", background_dsigma_wwj_production, "royalblue", (0, (5, 2))),
    (r"$ZZj$ $(\times 10^2)$", background_dsigma_zzj_production, "crimson", (0, (1, 1))),
    (r"$WZj$", background_dsigma_wzj_production, "slateblue", (0, (3, 2, 1, 2)))
]


for label, dsigma_dict, color, style in background_styles:
    m_w_bins, dsigma = dsigma_dict["m_w_leptonic_bins"]
    if np.sum(dsigma) > 0:
        normalized_dsigma = dsigma / np.sum(dsigma)
        plt.step(m_w_bins, normalized_dsigma, where="mid", alpha=0.7,
                 label=label, color=color, linestyle=style, linewidth=3)

# ✅ Axis Labels and Title
plt.xlabel(r"$M_W^{\ell \nu_\ell} \ \mathrm{[GeV]}$")
plt.ylabel("Normalized Entries")
plt.title(r"Delphes simulation : $e^- p \to e^- W^+ W^- p \to e^- j j \ell \nu_{\ell} p$ : LHeC@1.2 TeV", fontsize=20)
plt.ylim(0.0, 0.5)

# ✅ Legend and Layout
plt.legend()
plt.grid(True, linestyle="--", alpha=0.6)
plt.tight_layout()

# ✅ Save and Show
plt.savefig("/home/hamzeh-khanpour/Documents/GitHub/LHeC_Fast_Simulation/Pythia8_Delphes_semi_leptonic_allsignal_bkgs/normalized_m_w_leptonic_allFMsignal_allbkgs_all_decay_mode.pdf", dpi=600)
plt.show()






# ===================================================
# ✅ Normalized \eta_{leading \; jet} Distribution (Semi-Leptonic)
# ===================================================

plt.figure(figsize=(11, 12))
plt.subplots_adjust(left=0.15, right=0.95, bottom=0.12, top=0.95)

# ✅ Plot Normalized Signal Distributions
for signal_name, dsigma_data in signal_dsigma.items():
    centrality_bins, dsigma = dsigma_data["centrality_bins"]
    if np.sum(dsigma) > 0:
        normalized_dsigma = dsigma / np.sum(dsigma)
        plt.step(centrality_bins, normalized_dsigma, where="mid", alpha=0.7,
                 label=f"Signal ($W^+W^-$) [{signal_name}]",
                 color=signal_colors.get(signal_name, "black"), linewidth=3)

# ✅ Plot Normalized Backgrounds
background_styles = [
    (r"SM ($\gamma\gamma \to W^+W^-$)", background_dsigma_aa_ww, "blue", "-"),
    (r"$\gamma\gamma \to t\bar{t}$ $(\times 10^2)$", background_dsigma_aa_ttbar, "teal", "--"),
    (r"$\gamma\gamma \to \tau^+\tau^-$", background_dsigma_aa_tautau, "orange", "-"),
    (r"$\gamma\gamma \to \mu^+\mu^-$", background_dsigma_aa_mumu, "brown", (0, (3, 1, 1, 1))),
    (r"Inclusive $t\bar{t}$", background_dsigma_inclusive_ttbar, "orchid", (0, (3, 1, 1, 1))),
    (r"Single Top", background_dsigma_single_top, "teal", "--"),
    (r"W Production", background_dsigma_w_production, "darkgreen", "-."),
    (r"Z Production", background_dsigma_z_production, "darkred", ":"),
    (r"$WWj$", background_dsigma_wwj_production, "royalblue", (0, (5, 2))),
    (r"$ZZj$ $(\times 10^2)$", background_dsigma_zzj_production, "crimson", (0, (1, 1))),
    (r"$WZj$", background_dsigma_wzj_production, "slateblue", (0, (3, 2, 1, 2)))
]


for label, dsigma_dict, color, style in background_styles:
    centrality_bins, dsigma = dsigma_dict["centrality_bins"]
    if np.sum(dsigma) > 0:
        normalized_dsigma = dsigma / np.sum(dsigma)
        plt.step(centrality_bins, normalized_dsigma, where="mid", alpha=0.7,
                 label=label, color=color, linestyle=style, linewidth=3)

# ✅ Axis Labels and Title
plt.xlabel(r"$\eta_{leading \; jet}$")
plt.ylabel("Normalized Entries")
plt.title(r"Delphes simulation : $e^- p \to e^- W^+ W^- p \to e^- j j \ell \nu_{\ell} p$ : LHeC@1.2 TeV", fontsize=20)
plt.ylim(0.0, 0.2)

# ✅ Legend and Layout
plt.legend()
plt.grid(True, linestyle="--", alpha=0.6)
plt.tight_layout()

# ✅ Save and Show
plt.savefig("/home/hamzeh-khanpour/Documents/GitHub/LHeC_Fast_Simulation/Pythia8_Delphes_semi_leptonic_allsignal_bkgs/normalized_eta_leading_jet_allFMsignal_allbkgs_all_decay_mode.pdf", dpi=600)
plt.show()






# ===================================================
# ✅ Normalized eta_{subleading \; jet Distribution (Semi-Leptonic)
# ===================================================



plt.figure(figsize=(11, 12))
plt.subplots_adjust(left=0.15, right=0.95, bottom=0.12, top=0.95)

# ✅ Plot Normalized Signal Distributions
for signal_name, dsigma_data in signal_dsigma.items():
    exp_centrality_bins, dsigma = dsigma_data["exp_centrality_bins"]
    if np.sum(dsigma) > 0:
        normalized_dsigma = dsigma / np.sum(dsigma)
        plt.step(exp_centrality_bins, normalized_dsigma, where="mid", alpha=0.7,
                 label=f"Signal ($W^+W^-$) [{signal_name}]",
                 color=signal_colors.get(signal_name, "black"), linewidth=3)

# ✅ Plot Normalized Backgrounds
background_styles = [
    (r"SM ($\gamma\gamma \to W^+W^-$)", background_dsigma_aa_ww, "blue", "-"),
    (r"$\gamma\gamma \to t\bar{t}$ $(\times 10^2)$", background_dsigma_aa_ttbar, "teal", "--"),
    (r"$\gamma\gamma \to \tau^+\tau^-$", background_dsigma_aa_tautau, "orange", "-"),
    (r"$\gamma\gamma \to \mu^+\mu^-$", background_dsigma_aa_mumu, "brown", (0, (3, 1, 1, 1))),
    (r"Inclusive $t\bar{t}$", background_dsigma_inclusive_ttbar, "orchid", (0, (3, 1, 1, 1))),
    (r"Single Top", background_dsigma_single_top, "teal", "--"),
    (r"W Production", background_dsigma_w_production, "darkgreen", "-."),
    (r"Z Production", background_dsigma_z_production, "darkred", ":"),
    (r"$WWj$", background_dsigma_wwj_production, "royalblue", (0, (5, 2))),         # custom dashed
    (r"$ZZj$ $(\times 10^2)$", background_dsigma_zzj_production, "crimson", (0, (1, 1))),           # dense dots
    (r"$WZj$", background_dsigma_wzj_production, "slateblue", (0, (3, 2, 1, 2)))    # dash-dot-dash
]


for label, dsigma_dict, color, style in background_styles:
    exp_centrality_bins, dsigma = dsigma_dict["exp_centrality_bins"]
    if np.sum(dsigma) > 0:
        normalized_dsigma = dsigma / np.sum(dsigma)
        plt.step(exp_centrality_bins, normalized_dsigma, where="mid", alpha=0.7,
                 label=label, color=color, linestyle=style, linewidth=3)

# ✅ Axis Labels and Title
plt.xlabel(r"$\eta_{subleading \; jet}$")
plt.ylabel("Normalized Entries")
plt.title(r"Delphes simulation : $e^- p \to e^- W^+ W^- p \to e^- j j \ell \nu_{\ell} p$ : LHeC@1.2 TeV", fontsize=20)
plt.ylim(0.0, 0.2)

# ✅ Legend and Layout
plt.legend()
plt.grid(True, linestyle="--", alpha=0.6)
plt.tight_layout()

# ✅ Save and Show
plt.savefig("/home/hamzeh-khanpour/Documents/GitHub/LHeC_Fast_Simulation/Pythia8_Delphes_semi_leptonic_allsignal_bkgs/normalized_eta_subleading_jet_allFMsignal_allbkgs_all_decay_mode.pdf", dpi=600)
plt.show()





# ===================================================
# ✅ Normalized Jet Centrality Distribution (Semi-Leptonic)
# ===================================================

plt.figure(figsize=(11, 12))
plt.subplots_adjust(left=0.15, right=0.95, bottom=0.12, top=0.95)

# ✅ Plot Normalized Signal Distributions
for signal_name, dsigma_data in signal_dsigma.items():
    jet_centrality_bins, dsigma = dsigma_data["jet_centrality_bins"]
    if np.sum(dsigma) > 0:
        normalized_dsigma = dsigma / np.sum(dsigma)
        plt.step(jet_centrality_bins, normalized_dsigma, where="mid", alpha=0.7,
                 label=f"Signal ($W^+W^-$) [{signal_name}]",
                 color=signal_colors.get(signal_name, "black"), linewidth=3)

# ✅ Plot Normalized Backgrounds
background_styles = [
    (r"SM ($\gamma\gamma \to W^+W^-$)", background_dsigma_aa_ww, "blue", "-"),
    (r"$\gamma\gamma \to t\bar{t}$ $(\times 10^2)$", background_dsigma_aa_ttbar, "teal", "--"),
    (r"$\gamma\gamma \to \tau^+\tau^-$", background_dsigma_aa_tautau, "orange", "-"),
    (r"$\gamma\gamma \to \mu^+\mu^-$", background_dsigma_aa_mumu, "brown", (0, (3, 1, 1, 1))),
    (r"Inclusive $t\bar{t}$", background_dsigma_inclusive_ttbar, "orchid", (0, (3, 1, 1, 1))),
    (r"Single Top", background_dsigma_single_top, "teal", "--"),
    (r"W Production", background_dsigma_w_production, "darkgreen", "-."),
    (r"Z Production", background_dsigma_z_production, "darkred", ":"),
    (r"$WWj$", background_dsigma_wwj_production, "royalblue", (0, (5, 2))),         # new
    (r"$ZZj$ $(\times 10^2)$", background_dsigma_zzj_production, "crimson", (0, (1, 1))),           # new
    (r"$WZj$", background_dsigma_wzj_production, "slateblue", (0, (3, 2, 1, 2)))    # new
]


for label, dsigma_dict, color, style in background_styles:
    jet_centrality_bins, dsigma = dsigma_dict["jet_centrality_bins"]
    if np.sum(dsigma) > 0:
        normalized_dsigma = dsigma / np.sum(dsigma)
        plt.step(jet_centrality_bins, normalized_dsigma, where="mid", alpha=0.7,
                 label=label, color=color, linestyle=style, linewidth=3)

# ✅ Axis Labels and Title
plt.xlabel(r"$C_{\mathrm{jets}}$")
plt.ylabel("Normalized Entries")
plt.title(r"Delphes simulation : $e^- p \to e^- W^+ W^- p \to e^- j j \ell \nu_{\ell} p$ : LHeC@1.2 TeV", fontsize=20)
plt.ylim(0.0, 0.20)

# ✅ Legend and Layout
plt.legend()
plt.grid(True, linestyle="--", alpha=0.6)
# ✅ Add Jet Centrality Formula to Plot
plt.text(0.2, 0.15, r"$C_{\mathrm{jets}} = \frac{|\eta_{\mathrm{jet1}} + \eta_{\mathrm{jet2}}|}{2}$", fontsize=25, color="black")
plt.tight_layout()

# ✅ Save and Show
plt.savefig("/home/hamzeh-khanpour/Documents/GitHub/LHeC_Fast_Simulation/Pythia8_Delphes_semi_leptonic_allsignal_bkgs/normalized_jet_centrality_allFMsignal_allbkgs_all_decay_mode.pdf", dpi=600)
plt.show()







#=========================================================================
#=========================================================================

# Open the output ROOT file (use "UPDATE" if you want to append)
#output_file = ROOT.TFile("/home/hamzeh-khanpour/Documents/GitHub/LHeC_Fast_Simulation/Pythia8_Delphes_semi_leptonic_allsignal_bkgs/output_histograms.root", "RECREATE")

try:
    # ✅ Fix 1: Rename signal names to be ROOT-compatible (remove LaTeX)
    signal_dirs = {}
    for signal_name in signal_files.keys():
        clean_signal_name = signal_name.replace("$", "").replace("{", "").replace("}", "").replace("\\", "").replace(" ", "_")
        signal_dirs[signal_name] = output_file.mkdir(f"signal_{clean_signal_name}")

    # ✅ Fix 2: Create Background Directory BEFORE Writing Any Data
    background_dir = output_file.mkdir("SM_background")

    # ✅ Fix 3: Ensure `background_histograms` is properly defined
    background_histograms = {
        "hist_lepton": hist_lepton_background,
        "hist_leading_jet": hist_leading_jet_background,
        "hist_lepton_eta": hist_lepton_eta_background,
        "hist_delta_r": hist_delta_r_background,
        "hist_missing_et": hist_missing_et_background,
        "hist_subleading_jet_eta": hist_subleading_jet_eta_background,
        "hist_leading_jet_eta": hist_leading_jet_eta_background,
        "hist_jet_centrality": hist_jet_centrality_background,
        "hist_delta_eta_jj": hist_delta_eta_jj_background,
        "hist_m_w_hadronic": hist_m_w_hadronic_background,
        "hist_m_w_leptonic": hist_m_w_leptonic_background
    }

    # ✅ Fix 4: Save signal histograms (DO NOT redefine signal_histograms)
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

    print("✅ Histograms successfully saved to output_histograms.root with separate branches for each signal and 'SM_background'.")

except Exception as e:
    print(f"❌ Error while saving histograms: {e}")

finally:
    # ✅ Fix 6: Always close the ROOT file properly
    output_file.Close()
    print("📁 ROOT file closed successfully.")



#=========================================================================
#=========================================================================



