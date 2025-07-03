

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
#    "$FM_{0} / \Lambda^4$": "aa_ww_semi_leptonic_NP_1_FM0_Delphes_Pythia.root",
#    "$FM_{1} / \Lambda^4$": "aa_ww_semi_leptonic_NP_1_FM1_Delphes_Pythia.root",
    "FM2_Lambda4": "aa_ww_semi_leptonic_NP_1_FM1_Delphes_Pythia.root",
#    "$FM_{3} / \Lambda^4$": "aa_ww_semi_leptonic_NP_1_FM3_Delphes_Pythia.root",
}



aa_ww_background_file_path = "aa_ww_semi_leptonic_SM_NP_0_FMi_0_Delphes_Pythia.root"   # SM aa to WW


aa_ttbar_background_file_path  = "aa_ttbar_inclusive_decay_Delphes_Pythia.root"
#aa_tautau_background_file_path = "aa_tautau_Delphes_Pythia.root"
aa_tautau_background_file_path = "cepgen-gammagammatotautau-pt2p5-eta5_tauola.root"   # Laurent

aa_tautau_inel_background_file_path   = "cepgen-gammagammatotautau_sd_luxlike-pt2p5-eta5_tauola-pythia6.root"  # Laurent



inclusive_ttbar_background_file_path   = "ttbar_inclusive_decay.root"
single_top_background_file_path        = "single_top_inclusive_decay.root"


w_production_background_file_path      = "w_production_inclusive_decay_Delphes_Pythia.root"
z_production_background_file_path      = "z_production_inclusive_decay_Delphes_Pythia.root"


wwj_production_background_file_path      = "diboson_wwj_inclusive_decay_Delphes_Pythia.root"
zzj_production_background_file_path      = "diboson_zzj_inclusive_decay_Delphes_Pythia.root"
wzj_production_background_file_path      = "diboson_wzj_inclusive_decay_Delphes_Pythia.root"


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
    hist_pt_w_leptonic,
    hist_pt_w_hadronic,
    hist_delta_phi_lep_met,
    hist_mt_w_leptonic,
    hist_ht_total,
    hist_delta_phi_jj,
    hist_delta_phi_wl_wh,
    hist_delta_eta_wl_wh,
    hist_m_jj,
    hist_m_lvjj,
):



    # Open the ROOT file
    chain = ROOT.TChain("Delphes")
    chain.Add(file_path)

    # Create ExRootTreeReader object
    treeReader = ROOT.ExRootTreeReader(chain)
    numberOfEntries  =  treeReader.GetEntries()

    # Counters for efficiency calculation
    # ✅ Initialize selection counters
    total_events = numberOfEntries

    selected_events_pre = 0   # ✅ Fix: Initialize before use
    selected_events_pre_lepton = 0
    selected_events_pre_jets = 0
    selected_events_pre_eta_lepton = 0
    selected_events_pre_jet_centrality = 0

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
            if electron.PT > 1:
                leptons.append(lepton_vec)


        for i in range(branchMuon.GetEntries()):
            muon = branchMuon.At(i)
            lepton_vec = TLorentzVector()
            lepton_vec.SetPtEtaPhiM(muon.PT, muon.Eta, muon.Phi, 0.0)
            if muon.PT > 1:
                leptons.append(lepton_vec)



        if hist_missing_et is not None and branchMissingET.GetEntries() > 0:
            missing_et = branchMissingET.At(0)



        # Count the number of jets
        jets = []
        for i in range(branchJet.GetEntries()):
            jet = branchJet.At(i)
            jet_vec = TLorentzVector()
            jet_vec.SetPtEtaPhiM(jet.PT, jet.Eta, jet.Phi, jet.Mass)
            if jet.PT > 1:
                jets.append(jet_vec)



       # Apply selection criteria: exactly one lepton and exactly two jets

        if len(leptons) != 1 or len(jets) != 2:
#            continue
        # Count selected events
            selected_events_pre += 1 # ✅ Count events passing both lepton and jet requirements


        # ----------------------------------------
        # Step 1: Lepton requirement (1 lepton + pT > 10 GeV)

        if len(leptons) != 1:
            continue
        if leptons[0].Pt() <= 10:
            continue

           # ✅ MET cut  :   # Skip event if MET is too low
        if missing_et.MET < 10:
            continue

        selected_events_pre_lepton += 1



        # ----------------------------------------

        # Step 2: Exactly 2 jets with pT > 10 GeV
        if len(jets) != 2:
            continue
        if jets[0].Pt() <= 10 or jets[1].Pt() <= 10:
            continue
        selected_events_pre_jets += 1




        # Fill histogram with lepton pT
        hist_lepton_pt.Fill(leptons[0].Pt())
        hist_lepton_eta.Fill(leptons[0].Eta())   # ✅ Fill η histogram

        # Fill histogram with leading jet pT
        leading_jet = jets[0] if jets[0].Pt() > jets[1].Pt() else jets[1]
        hist_leading_jet_pt.Fill(leading_jet.Pt())


        # Fill histogram with missing_et.MET
        hist_missing_et.Fill(missing_et.MET)



        # Fill optional histograms  Fill histogram with leading jet pT
        if hist_delta_r is not None:
            delta_r = leptons[0].DeltaR(leading_jet)
            hist_delta_r.Fill(delta_r)


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
            delta_eta_jj = (jets[0].Eta() - jets[1].Eta())
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



            # ===============================================
            # ✅ New Observables for Semi-Leptonic WW Analysis
            # ===============================================

            if 'hist_pt_w_leptonic' in locals() and hist_pt_w_leptonic is not None:
                hist_pt_w_leptonic.Fill(w_leptonic.Pt())


            if 'hist_pt_w_hadronic' in locals() and hist_pt_w_hadronic is not None:
                hist_pt_w_hadronic.Fill(w_hadronic.Pt())


            if 'hist_delta_phi_lep_met' in locals() and hist_delta_phi_lep_met is not None:
                delta_phi = abs(leptons[0].Phi() - missing_et.Phi)
                if delta_phi > np.pi:
                    delta_phi = 2 * np.pi - delta_phi
                hist_delta_phi_lep_met.Fill(delta_phi)


            if 'hist_mt_w_leptonic' in locals() and hist_mt_w_leptonic is not None:
                mt_lep = np.sqrt(2 * leptons[0].Pt() * missing_et.MET * (1 - np.cos(delta_phi)))
                hist_mt_w_leptonic.Fill(mt_lep)


            if 'hist_ht_total' in locals() and hist_ht_total is not None:
                ht_total = leptons[0].Pt() + jets[0].Pt() + jets[1].Pt()
                hist_ht_total.Fill(ht_total)




            # =========================================
            # ✅ Extended Kinematic Observables (Optional)
            # =========================================

            if 'hist_delta_phi_jj' in locals() and hist_delta_phi_jj is not None:
                delta_phi_jj = (jets[0].Phi() - jets[1].Phi())
                if delta_phi_jj > np.pi:
                    delta_phi_jj = 2 * np.pi - delta_phi_jj
                hist_delta_phi_jj.Fill(delta_phi_jj)


            if 'hist_delta_phi_wl_wh' in locals() and hist_delta_phi_wl_wh is not None:
                delta_phi_wl_wh = (w_leptonic.Phi() - w_hadronic.Phi())
                if delta_phi_wl_wh > np.pi:
                    delta_phi_wl_wh = 2 * np.pi - delta_phi_wl_wh
                hist_delta_phi_wl_wh.Fill(delta_phi_wl_wh)


            if 'hist_delta_eta_wl_wh' in locals() and hist_delta_eta_wl_wh is not None:
                delta_eta_wl_wh = w_leptonic.Eta() - w_hadronic.Eta()
                hist_delta_eta_wl_wh.Fill(delta_eta_wl_wh)


            if 'hist_m_jj' in locals() and hist_m_jj is not None:
                hist_m_jj.Fill(w_hadronic.M())  # Already computed as dijet mass


            if 'hist_m_lvjj' in locals() and hist_m_lvjj is not None:
                lvjj = leptons[0] + neutrino_vec + jets[0] + jets[1]
                hist_m_lvjj.Fill(lvjj.M())




        #if abs(leptons[0].Eta()) > 0.5:
            #selected_events_pre_eta_lepton += 1

        #if jet_centrality > 2.0:
            #selected_events_pre_jet_centrality += 1

        ## **✅ Final Event Selection (W mass window)**
        #if 65 < w_leptonic.M() < 95 and 65 < w_hadronic.M() < 95:
           #selected_events_final += 1     # Count selected events


# ✅ Apply cuts in sequence: η_lepton → jet_centrality → W mass window
        if abs(leptons[0].Eta()) > -10000.0:    # 0.5:
            selected_events_pre_eta_lepton += 1

            if jet_centrality > -10000.0:   #   2.0:
                selected_events_pre_jet_centrality += 1

                if 65 < w_leptonic.M() < 95 and 65 < w_hadronic.M() < 95:
                    selected_events_final += 1  # ✅ Count event only if it passed all previous cuts



    # **✅ Selection Efficiency Calculations**
    efficiency_pre = selected_events_pre / total_events if total_events > 0 else 0

    efficiency_pre_lepton = selected_events_pre_lepton / total_events if total_events > 0 else 0
    efficiency_pre_jets = selected_events_pre_jets / total_events if total_events > 0 else 0


    efficiency_eta_lepton = selected_events_pre_eta_lepton / total_events if total_events > 0 else 0
    efficiency_jet_centrality = selected_events_pre_jet_centrality / total_events if total_events > 0 else 0


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
        hist_pt_w_leptonic,
        hist_pt_w_hadronic,
        hist_delta_phi_lep_met,
        hist_mt_w_leptonic,
        hist_ht_total,
        hist_delta_phi_jj,
        hist_delta_phi_wl_wh,
        hist_delta_eta_wl_wh,
        hist_m_jj,
        hist_m_lvjj,
        efficiency_pre_lepton,
        efficiency_pre_jets,
        efficiency_eta_lepton,
        efficiency_jet_centrality,
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
    "$FM_{0} / \Lambda^4$": 0.015395800000000001,    # pb
    "$FM_{1} / \Lambda^4$": 0.015587299999999998,    # pb
    "FM2_Lambda4": 0.022850099999999998,             # pb
    "$FM_{3} / \Lambda^4$": 0.017148100000000000     # pb
}


aa_ww_background_cross_section = 0.0149219 # aa_ww_semi_leptonic_SM_NP_0_FMi_0


aa_ttbar_background_cross_section  = 4.824774e-05 * 100.0  # pb  * 10^{+2}

#aa_tautau_background_cross_section = 14.46200  # pb
aa_tautau_background_cross_section   = 1.806765e-01  # pb  Laurent

aa_tautau_inel_background_cross_section = 1.156165e-01  # pb




inclusive_ttbar_background_cross_section   = 0.00817326           # pb
single_top_background_cross_section        = 1.36211000           # pb


w_production_background_cross_section      = 1.965201542           # pb
z_production_background_cross_section      = 0.159347434           # pb


wwj_production_background_cross_section      = 0.02031491612402401   # pb
zzj_production_background_cross_section      = 8.106588466651764e-05 * 100.0  # pb  * 10^{+2}
wzj_production_background_cross_section      = 0.0028587542003382592  # pb



num_bins = 50

# Existing ranges
pt_range_lepton         = (0, 400)     # Lepton pT
pt_range_jet            = (0, 400)     # Jet pT
eta_range               = (-4, 6)      # Pseudorapidity
delta_r_range           = (0, 9)       # ΔR
met_range               = (0, 300)     # MET
centrality_range        = (-3, 6)      # Centrality
exp_centrality_range    = (-3, 6)      # Exp centrality
jet_centrality_range    = (-1, 6)      # Jet centrality
delta_eta_jj_range      = (-6, 6)      # Δη(jj)

# W mass ranges
m_w_hadronic_range      = (1, 140)     # Hadronic W mass
m_w_leptonic_range      = (1, 140)     # Leptonic W mass

# New variable ranges
pt_range_w              = (0, 500)     # pT(W) for both leptonic and hadronic
delta_phi_range         = (0, np.pi)   # Δφ
mt_w_leptonic_range     = (0, 300)     # Transverse mass of W_lep
ht_total_range          = (0, 1000)    # Total HT

# Optional extensions
delta_phi_jj_range      = (0, np.pi)   # Δφ(jj)
delta_phi_wl_wh_range   = (0, np.pi)   # Δφ(W_lep, W_had)
delta_eta_wl_wh_range   = (-10, 10)    # Δη(W_lep, W_had)
m_jj_range              = (1, 140)     # m(jj) ~ m_W_had
m_lvjj_range            = (0, 1000)    # m(ℓνjj)

# Bin widths
bin_width_pt_lepton         = (pt_range_lepton[1] - pt_range_lepton[0]) / num_bins
bin_width_pt_jet            = (pt_range_jet[1] - pt_range_jet[0]) / num_bins
bin_width_eta               = (eta_range[1] - eta_range[0]) / num_bins
bin_width_delta_r           = (delta_r_range[1] - delta_r_range[0]) / num_bins
bin_width_met               = (met_range[1] - met_range[0]) / num_bins
bin_width_centrality        = (centrality_range[1] - centrality_range[0]) / num_bins
bin_width_exp_centrality    = (exp_centrality_range[1] - exp_centrality_range[0]) / num_bins
bin_width_jet_centrality    = (jet_centrality_range[1] - jet_centrality_range[0]) / num_bins
bin_width_delta_eta_jj      = (delta_eta_jj_range[1] - delta_eta_jj_range[0]) / num_bins
bin_width_m_w_hadronic      = (m_w_hadronic_range[1] - m_w_hadronic_range[0]) / num_bins
bin_width_m_w_leptonic      = (m_w_leptonic_range[1] - m_w_leptonic_range[0]) / num_bins

# Bin widths for new observables
bin_width_pt_w              = (pt_range_w[1] - pt_range_w[0]) / num_bins
bin_width_delta_phi         = (delta_phi_range[1] - delta_phi_range[0]) / num_bins
bin_width_mt_w_leptonic     = (mt_w_leptonic_range[1] - mt_w_leptonic_range[0]) / num_bins
bin_width_ht_total          = (ht_total_range[1] - ht_total_range[0]) / num_bins
bin_width_delta_phi_jj      = (delta_phi_jj_range[1] - delta_phi_jj_range[0]) / num_bins
bin_width_delta_phi_wl_wh   = (delta_phi_wl_wh_range[1] - delta_phi_wl_wh_range[0]) / num_bins
bin_width_delta_eta_wl_wh   = (delta_eta_wl_wh_range[1] - delta_eta_wl_wh_range[0]) / num_bins
bin_width_m_jj              = (m_jj_range[1] - m_jj_range[0]) / num_bins
bin_width_m_lvjj            = (m_lvjj_range[1] - m_lvjj_range[0]) / num_bins





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
        "hist_lepton_pt": ROOT.TH1F(f"hist_lepton_pt_{signal_name}", f"Lepton pT Distribution ({signal_name}); p_{{T}} [GeV]; Entries", num_bins, *pt_range_lepton),
        "hist_leading_jet_pt": ROOT.TH1F(f"hist_leading_jet_pt_{signal_name}", f"Leading Jet pT Distribution ({signal_name}); p_{{T}} [GeV]; Entries", num_bins, *pt_range_jet),
        "hist_lepton_eta": ROOT.TH1F(f"hist_lepton_eta_{signal_name}", f"Lepton Eta Distribution ({signal_name}); #eta; Entries", num_bins, *eta_range),
        "hist_delta_r": ROOT.TH1F(f"hist_delta_r_{signal_name}", f"Delta R Distribution ({signal_name}); ΔR; Entries", num_bins, *delta_r_range),
        "hist_missing_et": ROOT.TH1F(f"hist_missing_et_{signal_name}", f"Missing ET Distribution ({signal_name}); MET [GeV]; Entries", num_bins, *met_range),
        "hist_subleading_jet_eta": ROOT.TH1F(f"hist_subleading_jet_eta_{signal_name}", f"Centrality Distribution ({signal_name}); Centrality; Entries", num_bins, *centrality_range),
        "hist_leading_jet_eta": ROOT.TH1F(f"hist_leading_jet_eta_{signal_name}", f"Exponential Centrality Distribution ({signal_name}); exp(Centrality); Entries", num_bins, *exp_centrality_range),
        "hist_jet_centrality": ROOT.TH1F(f"hist_jet_centrality_{signal_name}", f"Jet Centrality Distribution ({signal_name}); Jet Centrality; Entries", num_bins, *jet_centrality_range),
        "hist_delta_eta_jj": ROOT.TH1F(f"hist_delta_eta_jj_{signal_name}", f"Delta Eta Between Jets ({signal_name}); Δηjj; Entries", num_bins, *delta_eta_jj_range),
        "hist_m_w_leptonic": ROOT.TH1F(f"hist_m_w_leptonic_{signal_name}", f"Leptonic W Boson Mass ({signal_name}); m_{{W}}^{{leptonic}} [GeV]; Entries", num_bins, *m_w_leptonic_range),
        "hist_m_w_hadronic": ROOT.TH1F(f"hist_m_w_hadronic_{signal_name}", f"Hadronic W Boson Mass ({signal_name}); m_{{W}}^{{hadronic}} [GeV]; Entries", num_bins, *m_w_hadronic_range),

        # ✅ New histograms for extended observables
        "hist_pt_w_leptonic": ROOT.TH1F(f"hist_pt_w_leptonic_{signal_name}", f"pT of Leptonic W ({signal_name}); p_{{T}} [GeV]; Entries", num_bins, *pt_range_w),
        "hist_pt_w_hadronic": ROOT.TH1F(f"hist_pt_w_hadronic_{signal_name}", f"pT of Hadronic W ({signal_name}); p_{{T}} [GeV]; Entries", num_bins, *pt_range_w),
        "hist_delta_phi_lep_met": ROOT.TH1F(f"hist_delta_phi_lep_met_{signal_name}", f"#Delta#phi(Lepton, MET) ({signal_name}); #Delta#phi; Entries", num_bins, *delta_phi_range),
        "hist_mt_w_leptonic": ROOT.TH1F(f"hist_mt_w_leptonic_{signal_name}", f"Transverse Mass of Leptonic W ({signal_name}); M_T [GeV]; Entries", num_bins, *mt_w_leptonic_range),
        "hist_ht_total": ROOT.TH1F(f"hist_ht_total_{signal_name}", f"Total H_T ({signal_name}); H_T [GeV]; Entries", num_bins, *ht_total_range),
        "hist_delta_phi_jj": ROOT.TH1F(f"hist_delta_phi_jj_{signal_name}", f"#Delta#phi(j_1, j_2) ({signal_name}); #Delta#phi_{{jj}}; Entries", num_bins, *delta_phi_jj_range),
        "hist_delta_phi_wl_wh": ROOT.TH1F(f"hist_delta_phi_wl_wh_{signal_name}", f"#Delta#phi(W^{{lep}}, W^{{had}}) ({signal_name}); #Delta#phi_{{WW}}; Entries", num_bins, *delta_phi_wl_wh_range),
        "hist_delta_eta_wl_wh": ROOT.TH1F(f"hist_delta_eta_wl_wh_{signal_name}", f"#Delta#eta(W^{{lep}}, W^{{had}}) ({signal_name}); #Delta#eta_{{WW}}; Entries", num_bins, *delta_eta_wl_wh_range),
        "hist_m_jj": ROOT.TH1F(f"hist_m_jj_{signal_name}", f"Dijet Mass m_{{jj}} ({signal_name}); m_{{jj}} [GeV]; Entries", num_bins, *m_jj_range),
        "hist_m_lvjj": ROOT.TH1F(f"hist_m_lvjj_{signal_name}", f"Invariant Mass m(ℓ#nu jj) ({signal_name}); m_{{ℓ#nu jj}} [GeV]; Entries", num_bins, *m_lvjj_range),
    }






# === AA_WW Background ===
hist_lepton_pt_aa_ww = ROOT.TH1F("hist_lepton_pt_aa_ww", "Lepton pT (aa_ww); p_{T} [GeV]; Entries", num_bins, *pt_range_lepton)
hist_leading_jet_pt_aa_ww = ROOT.TH1F("hist_leading_jet_pt_aa_ww", "Leading Jet pT (aa_ww); p_{T} [GeV]; Entries", num_bins, *pt_range_jet)
hist_lepton_eta_aa_ww = ROOT.TH1F("hist_lepton_eta_aa_ww", "Lepton Eta (aa_ww); #eta; Entries", num_bins, *eta_range)
hist_delta_r_aa_ww = ROOT.TH1F("hist_delta_r_aa_ww", "Delta R (aa_ww); ΔR; Entries", num_bins, *delta_r_range)
hist_missing_et_aa_ww = ROOT.TH1F("hist_missing_et_aa_ww", "MET (aa_ww); MET [GeV]; Entries", num_bins, *met_range)
hist_subleading_jet_eta_aa_ww = ROOT.TH1F("hist_subleading_jet_eta_aa_ww", "Centrality (aa_ww); Centrality; Entries", num_bins, *centrality_range)
hist_leading_jet_eta_aa_ww = ROOT.TH1F("hist_leading_jet_eta_aa_ww", "Exp Centrality (aa_ww); exp(Centrality); Entries", num_bins, *exp_centrality_range)
hist_jet_centrality_aa_ww = ROOT.TH1F("hist_jet_centrality_aa_ww", "Jet Centrality (aa_ww); Jet Centrality; Entries", num_bins, *jet_centrality_range)
hist_delta_eta_jj_aa_ww = ROOT.TH1F("hist_delta_eta_jj_aa_ww", "Delta Eta jj (aa_ww); Δηjj; Entries", num_bins, *delta_eta_jj_range)
hist_m_w_hadronic_aa_ww = ROOT.TH1F("hist_m_w_hadronic_aa_ww", "Hadronic W Mass (aa_ww); m_{W}^{hadronic} [GeV]; Entries", num_bins, *m_w_hadronic_range)
hist_m_w_leptonic_aa_ww = ROOT.TH1F("hist_m_w_leptonic_aa_ww", "Leptonic W Mass (aa_ww); m_{W}^{leptonic} [GeV]; Entries", num_bins, *m_w_leptonic_range)
# 10 new observables
hist_pt_w_leptonic_aa_ww = ROOT.TH1F("hist_pt_w_leptonic_aa_ww", "pT(W^leptonic) (aa_ww); p_{T} [GeV]; Entries", num_bins, 0, 400)
hist_pt_w_hadronic_aa_ww = ROOT.TH1F("hist_pt_w_hadronic_aa_ww", "pT(W^hadronic) (aa_ww); p_{T} [GeV]; Entries", num_bins, 0, 400)
hist_delta_phi_lep_met_aa_ww = ROOT.TH1F("hist_delta_phi_lep_met_aa_ww", "Delta Phi (lep, MET) (aa_ww); Δφ(ℓ, MET); Entries", num_bins, 0, np.pi)
hist_mt_w_leptonic_aa_ww = ROOT.TH1F("hist_mt_w_leptonic_aa_ww", "MT(W^leptonic) (aa_ww); M_{T} [GeV]; Entries", num_bins, 0, 400)
hist_ht_total_aa_ww = ROOT.TH1F("hist_ht_total_aa_ww", "HT Total (aa_ww); HT [GeV]; Entries", num_bins, 0, 1000)
hist_delta_phi_jj_aa_ww = ROOT.TH1F("hist_delta_phi_jj_aa_ww", "Δφ(j1, j2) (aa_ww); Δφ; Entries", num_bins, 0, np.pi)
hist_delta_phi_wl_wh_aa_ww = ROOT.TH1F("hist_delta_phi_wl_wh_aa_ww", "Δφ(Wlep, Whad) (aa_ww); Δφ; Entries", num_bins, 0, np.pi)
hist_delta_eta_wl_wh_aa_ww = ROOT.TH1F("hist_delta_eta_wl_wh_aa_ww", "Δη(Wlep, Whad) (aa_ww); Δη; Entries", num_bins, -6, 6)
hist_m_jj_aa_ww = ROOT.TH1F("hist_m_jj_aa_ww", "Dijet Invariant Mass (aa_ww); m_{jj} [GeV]; Entries", num_bins, 0, 400)
hist_m_lvjj_aa_ww = ROOT.TH1F("hist_m_lvjj_aa_ww", "m(lνjj) (aa_ww); m_{lνjj} [GeV]; Entries", num_bins, 0, 1000)







# === AA_TTBAR ===
hist_lepton_pt_aa_ttbar = ROOT.TH1F("hist_lepton_pt_aa_ttbar", "Lepton pT (aa_ttbar); p_{T} [GeV]; Entries", num_bins, *pt_range_lepton)
hist_leading_jet_pt_aa_ttbar = ROOT.TH1F("hist_leading_jet_pt_aa_ttbar", "Leading Jet pT (aa_ttbar); p_{T} [GeV]; Entries", num_bins, *pt_range_jet)
hist_lepton_eta_aa_ttbar = ROOT.TH1F("hist_lepton_eta_aa_ttbar", "Lepton Eta (aa_ttbar); #eta; Entries", num_bins, *eta_range)
hist_delta_r_aa_ttbar = ROOT.TH1F("hist_delta_r_aa_ttbar", "Delta R (aa_ttbar); ΔR; Entries", num_bins, *delta_r_range)
hist_missing_et_aa_ttbar = ROOT.TH1F("hist_missing_et_aa_ttbar", "MET (aa_ttbar); MET [GeV]; Entries", num_bins, *met_range)
hist_subleading_jet_eta_aa_ttbar = ROOT.TH1F("hist_subleading_jet_eta_aa_ttbar", "Centrality (aa_ttbar); Centrality; Entries", num_bins, *centrality_range)
hist_leading_jet_eta_aa_ttbar = ROOT.TH1F("hist_leading_jet_eta_aa_ttbar", "Exp Centrality (aa_ttbar); exp(Centrality); Entries", num_bins, *exp_centrality_range)
hist_jet_centrality_aa_ttbar = ROOT.TH1F("hist_jet_centrality_aa_ttbar", "Jet Centrality (aa_ttbar); Jet Centrality; Entries", num_bins, *jet_centrality_range)
hist_delta_eta_jj_aa_ttbar = ROOT.TH1F("hist_delta_eta_jj_aa_ttbar", "Delta Eta jj (aa_ttbar); Δηjj; Entries", num_bins, *delta_eta_jj_range)
hist_m_w_hadronic_aa_ttbar = ROOT.TH1F("hist_m_w_hadronic_aa_ttbar", "Hadronic W Mass (aa_ttbar); m_{W}^{hadronic} [GeV]; Entries", num_bins, *m_w_hadronic_range)
hist_m_w_leptonic_aa_ttbar = ROOT.TH1F("hist_m_w_leptonic_aa_ttbar", "Leptonic W Mass (aa_ttbar); m_{W}^{leptonic} [GeV]; Entries", num_bins, *m_w_leptonic_range)
# 10 new observables
hist_pt_w_leptonic_aa_ttbar = ROOT.TH1F("hist_pt_w_leptonic_aa_ttbar", "pT(W^leptonic) (aa_ttbar); p_{T} [GeV]; Entries", num_bins, 0, 400)
hist_pt_w_hadronic_aa_ttbar = ROOT.TH1F("hist_pt_w_hadronic_aa_ttbar", "pT(W^hadronic) (aa_ttbar); p_{T} [GeV]; Entries", num_bins, 0, 400)
hist_delta_phi_lep_met_aa_ttbar = ROOT.TH1F("hist_delta_phi_lep_met_aa_ttbar", "Delta Phi (lep, MET) (aa_ttbar); Δφ(ℓ, MET); Entries", num_bins, 0, np.pi)
hist_mt_w_leptonic_aa_ttbar = ROOT.TH1F("hist_mt_w_leptonic_aa_ttbar", "MT(W^leptonic) (aa_ttbar); M_{T} [GeV]; Entries", num_bins, 0, 400)
hist_ht_total_aa_ttbar = ROOT.TH1F("hist_ht_total_aa_ttbar", "HT Total (aa_ttbar); HT [GeV]; Entries", num_bins, 0, 1000)
hist_delta_phi_jj_aa_ttbar = ROOT.TH1F("hist_delta_phi_jj_aa_ttbar", "Δφ(j1, j2) (aa_ttbar); Δφ; Entries", num_bins, 0, np.pi)
hist_delta_phi_wl_wh_aa_ttbar = ROOT.TH1F("hist_delta_phi_wl_wh_aa_ttbar", "Δφ(Wlep, Whad) (aa_ttbar); Δφ; Entries", num_bins, 0, np.pi)
hist_delta_eta_wl_wh_aa_ttbar = ROOT.TH1F("hist_delta_eta_wl_wh_aa_ttbar", "Δη(Wlep, Whad) (aa_ttbar); Δη; Entries", num_bins, -6, 6)
hist_m_jj_aa_ttbar = ROOT.TH1F("hist_m_jj_aa_ttbar", "Dijet Invariant Mass (aa_ttbar); m_{jj} [GeV]; Entries", num_bins, 0, 400)
hist_m_lvjj_aa_ttbar = ROOT.TH1F("hist_m_lvjj_aa_ttbar", "m(lνjj) (aa_ttbar); m_{lνjj} [GeV]; Entries", num_bins, 0, 1000)










# === AA_TAUTAU ===
hist_lepton_pt_aa_tautau = ROOT.TH1F("hist_lepton_pt_aa_tautau", "Lepton pT (aa_tautau); p_{T} [GeV]; Entries", num_bins, *pt_range_lepton)
hist_leading_jet_pt_aa_tautau = ROOT.TH1F("hist_leading_jet_pt_aa_tautau", "Leading Jet pT (aa_tautau); p_{T} [GeV]; Entries", num_bins, *pt_range_jet)
hist_lepton_eta_aa_tautau = ROOT.TH1F("hist_lepton_eta_aa_tautau", "Lepton Eta (aa_tautau); #eta; Entries", num_bins, *eta_range)
hist_delta_r_aa_tautau = ROOT.TH1F("hist_delta_r_aa_tautau", "Delta R (aa_tautau); ΔR; Entries", num_bins, *delta_r_range)
hist_missing_et_aa_tautau = ROOT.TH1F("hist_missing_et_aa_tautau", "MET (aa_tautau); MET [GeV]; Entries", num_bins, *met_range)
hist_subleading_jet_eta_aa_tautau = ROOT.TH1F("hist_subleading_jet_eta_aa_tautau", "Centrality (aa_tautau); Centrality; Entries", num_bins, *centrality_range)
hist_leading_jet_eta_aa_tautau = ROOT.TH1F("hist_leading_jet_eta_aa_tautau", "Exp Centrality (aa_tautau); exp(Centrality); Entries", num_bins, *exp_centrality_range)
hist_jet_centrality_aa_tautau = ROOT.TH1F("hist_jet_centrality_aa_tautau", "Jet Centrality (aa_tautau); Jet Centrality; Entries", num_bins, *jet_centrality_range)
hist_delta_eta_jj_aa_tautau = ROOT.TH1F("hist_delta_eta_jj_aa_tautau", "Delta Eta jj (aa_tautau); Δηjj; Entries", num_bins, *delta_eta_jj_range)
hist_m_w_hadronic_aa_tautau = ROOT.TH1F("hist_m_w_hadronic_aa_tautau", "Hadronic W Mass (aa_tautau); m_{W}^{hadronic} [GeV]; Entries", num_bins, *m_w_hadronic_range)
hist_m_w_leptonic_aa_tautau = ROOT.TH1F("hist_m_w_leptonic_aa_tautau", "Leptonic W Mass (aa_tautau); m_{W}^{leptonic} [GeV]; Entries", num_bins, *m_w_leptonic_range)
# 10 new observables
hist_pt_w_leptonic_aa_tautau = ROOT.TH1F("hist_pt_w_leptonic_aa_tautau", "pT(W^leptonic) (aa_tautau); p_{T} [GeV]; Entries", num_bins, 0, 400)
hist_pt_w_hadronic_aa_tautau = ROOT.TH1F("hist_pt_w_hadronic_aa_tautau", "pT(W^hadronic) (aa_tautau); p_{T} [GeV]; Entries", num_bins, 0, 400)
hist_delta_phi_lep_met_aa_tautau = ROOT.TH1F("hist_delta_phi_lep_met_aa_tautau", "Delta Phi (lep, MET) (aa_tautau); Δφ(ℓ, MET); Entries", num_bins, 0, np.pi)
hist_mt_w_leptonic_aa_tautau = ROOT.TH1F("hist_mt_w_leptonic_aa_tautau", "MT(W^leptonic) (aa_tautau); M_{T} [GeV]; Entries", num_bins, 0, 400)
hist_ht_total_aa_tautau = ROOT.TH1F("hist_ht_total_aa_tautau", "HT Total (aa_tautau); HT [GeV]; Entries", num_bins, 0, 1000)
hist_delta_phi_jj_aa_tautau = ROOT.TH1F("hist_delta_phi_jj_aa_tautau", "Δφ(j1, j2) (aa_tautau); Δφ; Entries", num_bins, 0, np.pi)
hist_delta_phi_wl_wh_aa_tautau = ROOT.TH1F("hist_delta_phi_wl_wh_aa_tautau", "Δφ(Wlep, Whad) (aa_tautau); Δφ; Entries", num_bins, 0, np.pi)
hist_delta_eta_wl_wh_aa_tautau = ROOT.TH1F("hist_delta_eta_wl_wh_aa_tautau", "Δη(Wlep, Whad) (aa_tautau); Δη; Entries", num_bins, -6, 6)
hist_m_jj_aa_tautau = ROOT.TH1F("hist_m_jj_aa_tautau", "Dijet Invariant Mass (aa_tautau); m_{jj} [GeV]; Entries", num_bins, 0, 400)
hist_m_lvjj_aa_tautau = ROOT.TH1F("hist_m_lvjj_aa_tautau", "m(lνjj) (aa_tautau); m_{lνjj} [GeV]; Entries", num_bins, 0, 1000)







# === AA_tautau_inel ===
hist_lepton_pt_aa_tautau_inel = ROOT.TH1F("hist_lepton_pt_aa_tautau_inel", "Lepton pT (aa_tautau_inel); p_{T} [GeV]; Entries", num_bins, *pt_range_lepton)
hist_leading_jet_pt_aa_tautau_inel = ROOT.TH1F("hist_leading_jet_pt_aa_tautau_inel", "Leading Jet pT (aa_tautau_inel); p_{T} [GeV]; Entries", num_bins, *pt_range_jet)
hist_lepton_eta_aa_tautau_inel = ROOT.TH1F("hist_lepton_eta_aa_tautau_inel", "Lepton Eta (aa_tautau_inel); #eta; Entries", num_bins, *eta_range)
hist_delta_r_aa_tautau_inel = ROOT.TH1F("hist_delta_r_aa_tautau_inel", "Delta R (aa_tautau_inel); ΔR; Entries", num_bins, *delta_r_range)
hist_missing_et_aa_tautau_inel = ROOT.TH1F("hist_missing_et_aa_tautau_inel", "MET (aa_tautau_inel); MET [GeV]; Entries", num_bins, *met_range)
hist_subleading_jet_eta_aa_tautau_inel = ROOT.TH1F("hist_subleading_jet_eta_aa_tautau_inel", "Centrality (aa_tautau_inel); Centrality; Entries", num_bins, *centrality_range)
hist_leading_jet_eta_aa_tautau_inel = ROOT.TH1F("hist_leading_jet_eta_aa_tautau_inel", "Exp Centrality (aa_tautau_inel); exp(Centrality); Entries", num_bins, *exp_centrality_range)
hist_jet_centrality_aa_tautau_inel = ROOT.TH1F("hist_jet_centrality_aa_tautau_inel", "Jet Centrality (aa_tautau_inel); Jet Centrality; Entries", num_bins, *jet_centrality_range)
hist_delta_eta_jj_aa_tautau_inel = ROOT.TH1F("hist_delta_eta_jj_aa_tautau_inel", "Delta Eta jj (aa_tautau_inel); Δηjj; Entries", num_bins, *delta_eta_jj_range)
hist_m_w_hadronic_aa_tautau_inel = ROOT.TH1F("hist_m_w_hadronic_aa_tautau_inel", "Hadronic W Mass (aa_tautau_inel); m_{W}^{hadronic} [GeV]; Entries", num_bins, *m_w_hadronic_range)
hist_m_w_leptonic_aa_tautau_inel = ROOT.TH1F("hist_m_w_leptonic_aa_tautau_inel", "Leptonic W Mass (aa_tautau_inel); m_{W}^{leptonic} [GeV]; Entries", num_bins, *m_w_leptonic_range)
# 10 new observables
hist_pt_w_leptonic_aa_tautau_inel = ROOT.TH1F("hist_pt_w_leptonic_aa_tautau_inel", "pT(W^leptonic) (aa_tautau_inel); p_{T} [GeV]; Entries", num_bins, 0, 400)
hist_pt_w_hadronic_aa_tautau_inel = ROOT.TH1F("hist_pt_w_hadronic_aa_tautau_inel", "pT(W^hadronic) (aa_tautau_inel); p_{T} [GeV]; Entries", num_bins, 0, 400)
hist_delta_phi_lep_met_aa_tautau_inel = ROOT.TH1F("hist_delta_phi_lep_met_aa_tautau_inel", "Delta Phi (lep, MET) (aa_tautau_inel); Δφ(ℓ, MET); Entries", num_bins, 0, np.pi)
hist_mt_w_leptonic_aa_tautau_inel = ROOT.TH1F("hist_mt_w_leptonic_aa_tautau_inel", "MT(W^leptonic) (aa_tautau_inel); M_{T} [GeV]; Entries", num_bins, 0, 400)
hist_ht_total_aa_tautau_inel = ROOT.TH1F("hist_ht_total_aa_tautau_inel", "HT Total (aa_tautau_inel); HT [GeV]; Entries", num_bins, 0, 1000)
hist_delta_phi_jj_aa_tautau_inel = ROOT.TH1F("hist_delta_phi_jj_aa_tautau_inel", "Δφ(j1, j2) (aa_tautau_inel); Δφ; Entries", num_bins, 0, np.pi)
hist_delta_phi_wl_wh_aa_tautau_inel = ROOT.TH1F("hist_delta_phi_wl_wh_aa_tautau_inel", "Δφ(Wlep, Whad) (aa_tautau_inel); Δφ; Entries", num_bins, 0, np.pi)
hist_delta_eta_wl_wh_aa_tautau_inel = ROOT.TH1F("hist_delta_eta_wl_wh_aa_tautau_inel", "Δη(Wlep, Whad) (aa_tautau_inel); Δη; Entries", num_bins, -6, 6)
hist_m_jj_aa_tautau_inel = ROOT.TH1F("hist_m_jj_aa_tautau_inel", "Dijet Invariant Mass (aa_tautau_inel); m_{jj} [GeV]; Entries", num_bins, 0, 400)
hist_m_lvjj_aa_tautau_inel = ROOT.TH1F("hist_m_lvjj_aa_tautau_inel", "m(lνjj) (aa_tautau_inel); m_{lνjj} [GeV]; Entries", num_bins, 0, 1000)







# === Inclusive TTBAR ===
hist_lepton_pt_inclusive_ttbar             = ROOT.TH1F("hist_lepton_pt_inclusive_ttbar",             "Lepton pT (inclusive_ttbar); p_{T} [GeV]; Entries", num_bins, *pt_range_lepton)
hist_leading_jet_pt_inclusive_ttbar        = ROOT.TH1F("hist_leading_jet_pt_inclusive_ttbar",        "Leading Jet pT (inclusive_ttbar); p_{T} [GeV]; Entries", num_bins, *pt_range_jet)
hist_lepton_eta_inclusive_ttbar            = ROOT.TH1F("hist_lepton_eta_inclusive_ttbar",            "Lepton Eta (inclusive_ttbar); #eta; Entries", num_bins, *eta_range)
hist_delta_r_inclusive_ttbar               = ROOT.TH1F("hist_delta_r_inclusive_ttbar",               "Delta R (inclusive_ttbar); ΔR; Entries", num_bins, *delta_r_range)
hist_missing_et_inclusive_ttbar            = ROOT.TH1F("hist_missing_et_inclusive_ttbar",            "MET (inclusive_ttbar); MET [GeV]; Entries", num_bins, *met_range)
hist_subleading_jet_eta_inclusive_ttbar    = ROOT.TH1F("hist_subleading_jet_eta_inclusive_ttbar",    "Centrality (inclusive_ttbar); Centrality; Entries", num_bins, *centrality_range)
hist_leading_jet_eta_inclusive_ttbar       = ROOT.TH1F("hist_leading_jet_eta_inclusive_ttbar",       "Exp Centrality (inclusive_ttbar); exp(Centrality); Entries", num_bins, *exp_centrality_range)
hist_jet_centrality_inclusive_ttbar        = ROOT.TH1F("hist_jet_centrality_inclusive_ttbar",        "Jet Centrality (inclusive_ttbar); Jet Centrality; Entries", num_bins, *jet_centrality_range)
hist_delta_eta_jj_inclusive_ttbar          = ROOT.TH1F("hist_delta_eta_jj_inclusive_ttbar",          "Delta Eta jj (inclusive_ttbar); Δηjj; Entries", num_bins, *delta_eta_jj_range)
hist_m_w_hadronic_inclusive_ttbar          = ROOT.TH1F("hist_m_w_hadronic_inclusive_ttbar",          "Hadronic W Mass (inclusive_ttbar); m_{W}^{hadronic} [GeV]; Entries", num_bins, *m_w_hadronic_range)
hist_m_w_leptonic_inclusive_ttbar          = ROOT.TH1F("hist_m_w_leptonic_inclusive_ttbar",          "Leptonic W Mass (inclusive_ttbar); m_{W}^{leptonic} [GeV]; Entries", num_bins, *m_w_leptonic_range)
# === Additional Kinematic Observables ===
hist_pt_w_leptonic_inclusive_ttbar         = ROOT.TH1F("hist_pt_w_leptonic_inclusive_ttbar",         "pT(W^leptonic) (inclusive_ttbar); p_{T} [GeV]; Entries", num_bins, 0, 400)
hist_pt_w_hadronic_inclusive_ttbar         = ROOT.TH1F("hist_pt_w_hadronic_inclusive_ttbar",         "pT(W^hadronic) (inclusive_ttbar); p_{T} [GeV]; Entries", num_bins, 0, 400)
hist_delta_phi_lep_met_inclusive_ttbar     = ROOT.TH1F("hist_delta_phi_lep_met_inclusive_ttbar",     "Delta Phi (lep, MET) (inclusive_ttbar); Δφ(ℓ, MET); Entries", num_bins, 0, np.pi)
hist_mt_w_leptonic_inclusive_ttbar         = ROOT.TH1F("hist_mt_w_leptonic_inclusive_ttbar",         "MT(W^leptonic) (inclusive_ttbar); M_{T} [GeV]; Entries", num_bins, 0, 400)
hist_ht_total_inclusive_ttbar              = ROOT.TH1F("hist_ht_total_inclusive_ttbar",              "HT Total (inclusive_ttbar); HT [GeV]; Entries", num_bins, 0, 1000)
hist_delta_phi_jj_inclusive_ttbar          = ROOT.TH1F("hist_delta_phi_jj_inclusive_ttbar",          "Δφ(j1, j2) (inclusive_ttbar); Δφ; Entries", num_bins, 0, np.pi)
hist_delta_phi_wl_wh_inclusive_ttbar       = ROOT.TH1F("hist_delta_phi_wl_wh_inclusive_ttbar",       "Δφ(Wlep, Whad) (inclusive_ttbar); Δφ; Entries", num_bins, 0, np.pi)
hist_delta_eta_wl_wh_inclusive_ttbar       = ROOT.TH1F("hist_delta_eta_wl_wh_inclusive_ttbar",       "Δη(Wlep, Whad) (inclusive_ttbar); Δη; Entries", num_bins, -6, 6)
hist_m_jj_inclusive_ttbar                  = ROOT.TH1F("hist_m_jj_inclusive_ttbar",                  "Dijet Invariant Mass (inclusive_ttbar); m_{jj} [GeV]; Entries", num_bins, 0, 400)
hist_m_lvjj_inclusive_ttbar                = ROOT.TH1F("hist_m_lvjj_inclusive_ttbar",                "m(lνjj) (inclusive_ttbar); m_{lνjj} [GeV]; Entries", num_bins, 0, 1000)




# === SINGLE TOP ===
hist_lepton_pt_single_top             = ROOT.TH1F("hist_lepton_pt_single_top",             "Lepton pT (single_top); p_{T} [GeV]; Entries", num_bins, *pt_range_lepton)
hist_leading_jet_pt_single_top        = ROOT.TH1F("hist_leading_jet_pt_single_top",        "Leading Jet pT (single_top); p_{T} [GeV]; Entries", num_bins, *pt_range_jet)
hist_lepton_eta_single_top            = ROOT.TH1F("hist_lepton_eta_single_top",            "Lepton Eta (single_top); #eta; Entries", num_bins, *eta_range)
hist_delta_r_single_top               = ROOT.TH1F("hist_delta_r_single_top",               "Delta R (single_top); ΔR; Entries", num_bins, *delta_r_range)
hist_missing_et_single_top            = ROOT.TH1F("hist_missing_et_single_top",            "MET (single_top); MET [GeV]; Entries", num_bins, *met_range)
hist_subleading_jet_eta_single_top    = ROOT.TH1F("hist_subleading_jet_eta_single_top",    "Centrality (single_top); Centrality; Entries", num_bins, *centrality_range)
hist_leading_jet_eta_single_top       = ROOT.TH1F("hist_leading_jet_eta_single_top",       "Exp Centrality (single_top); exp(Centrality); Entries", num_bins, *exp_centrality_range)
hist_jet_centrality_single_top        = ROOT.TH1F("hist_jet_centrality_single_top",        "Jet Centrality (single_top); Jet Centrality; Entries", num_bins, *jet_centrality_range)
hist_delta_eta_jj_single_top          = ROOT.TH1F("hist_delta_eta_jj_single_top",          "Delta Eta jj (single_top); Δηjj; Entries", num_bins, *delta_eta_jj_range)
hist_m_w_hadronic_single_top          = ROOT.TH1F("hist_m_w_hadronic_single_top",          "Hadronic W Mass (single_top); m_{W}^{hadronic} [GeV]; Entries", num_bins, *m_w_hadronic_range)
hist_m_w_leptonic_single_top          = ROOT.TH1F("hist_m_w_leptonic_single_top",          "Leptonic W Mass (single_top); m_{W}^{leptonic} [GeV]; Entries", num_bins, *m_w_leptonic_range)
# === Additional Kinematic Observables ===
hist_pt_w_leptonic_single_top         = ROOT.TH1F("hist_pt_w_leptonic_single_top",         "pT(W^leptonic) (single_top); p_{T} [GeV]; Entries", num_bins, 0, 400)
hist_pt_w_hadronic_single_top         = ROOT.TH1F("hist_pt_w_hadronic_single_top",         "pT(W^hadronic) (single_top); p_{T} [GeV]; Entries", num_bins, 0, 400)
hist_delta_phi_lep_met_single_top     = ROOT.TH1F("hist_delta_phi_lep_met_single_top",     "Delta Phi (lep, MET) (single_top); Δφ(ℓ, MET); Entries", num_bins, 0, np.pi)
hist_mt_w_leptonic_single_top         = ROOT.TH1F("hist_mt_w_leptonic_single_top",         "MT(W^leptonic) (single_top); M_{T} [GeV]; Entries", num_bins, 0, 400)
hist_ht_total_single_top              = ROOT.TH1F("hist_ht_total_single_top",              "HT Total (single_top); HT [GeV]; Entries", num_bins, 0, 1000)
hist_delta_phi_jj_single_top          = ROOT.TH1F("hist_delta_phi_jj_single_top",          "Δφ(j1, j2) (single_top); Δφ; Entries", num_bins, 0, np.pi)
hist_delta_phi_wl_wh_single_top       = ROOT.TH1F("hist_delta_phi_wl_wh_single_top",       "Δφ(Wlep, Whad) (single_top); Δφ; Entries", num_bins, 0, np.pi)
hist_delta_eta_wl_wh_single_top       = ROOT.TH1F("hist_delta_eta_wl_wh_single_top",       "Δη(Wlep, Whad) (single_top); Δη; Entries", num_bins, -6, 6)
hist_m_jj_single_top                  = ROOT.TH1F("hist_m_jj_single_top",                  "Dijet Invariant Mass (single_top); m_{jj} [GeV]; Entries", num_bins, 0, 400)
hist_m_lvjj_single_top                = ROOT.TH1F("hist_m_lvjj_single_top",                "m(lνjj) (single_top); m_{lνjj} [GeV]; Entries", num_bins, 0, 1000)






# === W PRODUCTION ===
hist_lepton_pt_w_production             = ROOT.TH1F("hist_lepton_pt_w_production",             "Lepton pT (w_production); p_{T} [GeV]; Entries", num_bins, *pt_range_lepton)
hist_leading_jet_pt_w_production        = ROOT.TH1F("hist_leading_jet_pt_w_production",        "Leading Jet pT (w_production); p_{T} [GeV]; Entries", num_bins, *pt_range_jet)
hist_lepton_eta_w_production            = ROOT.TH1F("hist_lepton_eta_w_production",            "Lepton Eta (w_production); #eta; Entries", num_bins, *eta_range)
hist_delta_r_w_production               = ROOT.TH1F("hist_delta_r_w_production",               "Delta R (w_production); ΔR; Entries", num_bins, *delta_r_range)
hist_missing_et_w_production            = ROOT.TH1F("hist_missing_et_w_production",            "MET (w_production); MET [GeV]; Entries", num_bins, *met_range)
hist_subleading_jet_eta_w_production    = ROOT.TH1F("hist_subleading_jet_eta_w_production",    "Centrality (w_production); Centrality; Entries", num_bins, *centrality_range)
hist_leading_jet_eta_w_production       = ROOT.TH1F("hist_leading_jet_eta_w_production",       "Exp Centrality (w_production); exp(Centrality); Entries", num_bins, *exp_centrality_range)
hist_jet_centrality_w_production        = ROOT.TH1F("hist_jet_centrality_w_production",        "Jet Centrality (w_production); Jet Centrality; Entries", num_bins, *jet_centrality_range)
hist_delta_eta_jj_w_production          = ROOT.TH1F("hist_delta_eta_jj_w_production",          "Delta Eta jj (w_production); Δηjj; Entries", num_bins, *delta_eta_jj_range)
hist_m_w_hadronic_w_production          = ROOT.TH1F("hist_m_w_hadronic_w_production",          "Hadronic W Mass (w_production); m_{W}^{hadronic} [GeV]; Entries", num_bins, *m_w_hadronic_range)
hist_m_w_leptonic_w_production          = ROOT.TH1F("hist_m_w_leptonic_w_production",          "Leptonic W Mass (w_production); m_{W}^{leptonic} [GeV]; Entries", num_bins, *m_w_leptonic_range)
# === Additional Kinematic Observables ===
hist_pt_w_leptonic_w_production         = ROOT.TH1F("hist_pt_w_leptonic_w_production",         "pT(W^leptonic) (w_production); p_{T} [GeV]; Entries", num_bins, 0, 400)
hist_pt_w_hadronic_w_production         = ROOT.TH1F("hist_pt_w_hadronic_w_production",         "pT(W^hadronic) (w_production); p_{T} [GeV]; Entries", num_bins, 0, 400)
hist_delta_phi_lep_met_w_production     = ROOT.TH1F("hist_delta_phi_lep_met_w_production",     "Δφ(ℓ,MET) (w_production); Δφ; Entries", num_bins, 0, np.pi)
hist_mt_w_leptonic_w_production         = ROOT.TH1F("hist_mt_w_leptonic_w_production",         "MT(W^leptonic) (w_production); M_{T} [GeV]; Entries", num_bins, 0, 400)
hist_ht_total_w_production              = ROOT.TH1F("hist_ht_total_w_production",              "HT Total (w_production); H_{T} [GeV]; Entries", num_bins, 0, 1000)
hist_delta_phi_jj_w_production          = ROOT.TH1F("hist_delta_phi_jj_w_production",          "Δφ(j1, j2) (w_production); Δφ; Entries", num_bins, 0, np.pi)
hist_delta_phi_wl_wh_w_production       = ROOT.TH1F("hist_delta_phi_wl_wh_w_production",       "Δφ(Wlep, Whad) (w_production); Δφ; Entries", num_bins, 0, np.pi)
hist_delta_eta_wl_wh_w_production       = ROOT.TH1F("hist_delta_eta_wl_wh_w_production",       "Δη(Wlep, Whad) (w_production); Δη; Entries", num_bins, -6, 6)
hist_m_jj_w_production                  = ROOT.TH1F("hist_m_jj_w_production",                  "Dijet Invariant Mass (w_production); m_{jj} [GeV]; Entries", num_bins, 0, 400)
hist_m_lvjj_w_production                = ROOT.TH1F("hist_m_lvjj_w_production",                "m(lνjj) (w_production); m_{lνjj} [GeV]; Entries", num_bins, 0, 1000)








# === Z PRODUCTION ===
hist_lepton_pt_z_production             = ROOT.TH1F("hist_lepton_pt_z_production",             "Lepton pT (z_production); p_{T} [GeV]; Entries", num_bins, *pt_range_lepton)
hist_leading_jet_pt_z_production        = ROOT.TH1F("hist_leading_jet_pt_z_production",        "Leading Jet pT (z_production); p_{T} [GeV]; Entries", num_bins, *pt_range_jet)
hist_lepton_eta_z_production            = ROOT.TH1F("hist_lepton_eta_z_production",            "Lepton Eta (z_production); #eta; Entries", num_bins, *eta_range)
hist_delta_r_z_production               = ROOT.TH1F("hist_delta_r_z_production",               "Delta R (z_production); ΔR; Entries", num_bins, *delta_r_range)
hist_missing_et_z_production            = ROOT.TH1F("hist_missing_et_z_production",            "MET (z_production); MET [GeV]; Entries", num_bins, *met_range)
hist_subleading_jet_eta_z_production    = ROOT.TH1F("hist_subleading_jet_eta_z_production",    "Centrality (z_production); Centrality; Entries", num_bins, *centrality_range)
hist_leading_jet_eta_z_production       = ROOT.TH1F("hist_leading_jet_eta_z_production",       "Exp Centrality (z_production); exp(Centrality); Entries", num_bins, *exp_centrality_range)
hist_jet_centrality_z_production        = ROOT.TH1F("hist_jet_centrality_z_production",        "Jet Centrality (z_production); Jet Centrality; Entries", num_bins, *jet_centrality_range)
hist_delta_eta_jj_z_production          = ROOT.TH1F("hist_delta_eta_jj_z_production",          "Delta Eta jj (z_production); Δηjj; Entries", num_bins, *delta_eta_jj_range)
hist_m_w_hadronic_z_production          = ROOT.TH1F("hist_m_w_hadronic_z_production",          "Hadronic W Mass (z_production); m_{W}^{hadronic} [GeV]; Entries", num_bins, *m_w_hadronic_range)
hist_m_w_leptonic_z_production          = ROOT.TH1F("hist_m_w_leptonic_z_production",          "Leptonic W Mass (z_production); m_{W}^{leptonic} [GeV]; Entries", num_bins, *m_w_leptonic_range)
# === Additional Kinematic Observables ===
hist_pt_w_leptonic_z_production         = ROOT.TH1F("hist_pt_w_leptonic_z_production",         "pT(W^leptonic) (z_production); p_{T} [GeV]; Entries", num_bins, 0, 400)
hist_pt_w_hadronic_z_production         = ROOT.TH1F("hist_pt_w_hadronic_z_production",         "pT(W^hadronic) (z_production); p_{T} [GeV]; Entries", num_bins, 0, 400)
hist_delta_phi_lep_met_z_production     = ROOT.TH1F("hist_delta_phi_lep_met_z_production",     "Δφ(ℓ,MET) (z_production); Δφ; Entries", num_bins, 0, np.pi)
hist_mt_w_leptonic_z_production         = ROOT.TH1F("hist_mt_w_leptonic_z_production",         "MT(W^leptonic) (z_production); M_{T} [GeV]; Entries", num_bins, 0, 400)
hist_ht_total_z_production              = ROOT.TH1F("hist_ht_total_z_production",              "HT Total (z_production); H_{T} [GeV]; Entries", num_bins, 0, 1000)
hist_delta_phi_jj_z_production          = ROOT.TH1F("hist_delta_phi_jj_z_production",          "Δφ(j1, j2) (z_production); Δφ; Entries", num_bins, 0, np.pi)
hist_delta_phi_wl_wh_z_production       = ROOT.TH1F("hist_delta_phi_wl_wh_z_production",       "Δφ(Wlep, Whad) (z_production); Δφ; Entries", num_bins, 0, np.pi)
hist_delta_eta_wl_wh_z_production       = ROOT.TH1F("hist_delta_eta_wl_wh_z_production",       "Δη(Wlep, Whad) (z_production); Δη; Entries", num_bins, -6, 6)
hist_m_jj_z_production                  = ROOT.TH1F("hist_m_jj_z_production",                  "Dijet Invariant Mass (z_production); m_{jj} [GeV]; Entries", num_bins, 0, 400)
hist_m_lvjj_z_production                = ROOT.TH1F("hist_m_lvjj_z_production",                "m(lνjj) (z_production); m_{lνjj} [GeV]; Entries", num_bins, 0, 1000)







# === WWJ PRODUCTION ===
hist_lepton_pt_wwj_production             = ROOT.TH1F("hist_lepton_pt_wwj_production",             "Lepton pT (wwj_production); p_{T} [GeV]; Entries", num_bins, *pt_range_lepton)
hist_leading_jet_pt_wwj_production        = ROOT.TH1F("hist_leading_jet_pt_wwj_production",        "Leading Jet pT (wwj_production); p_{T} [GeV]; Entries", num_bins, *pt_range_jet)
hist_lepton_eta_wwj_production            = ROOT.TH1F("hist_lepton_eta_wwj_production",            "Lepton Eta (wwj_production); #eta; Entries", num_bins, *eta_range)
hist_delta_r_wwj_production               = ROOT.TH1F("hist_delta_r_wwj_production",               "Delta R (wwj_production); ΔR; Entries", num_bins, *delta_r_range)
hist_missing_et_wwj_production            = ROOT.TH1F("hist_missing_et_wwj_production",            "MET (wwj_production); MET [GeV]; Entries", num_bins, *met_range)
hist_subleading_jet_eta_wwj_production    = ROOT.TH1F("hist_subleading_jet_eta_wwj_production",    "Centrality (wwj_production); Centrality; Entries", num_bins, *centrality_range)
hist_leading_jet_eta_wwj_production       = ROOT.TH1F("hist_leading_jet_eta_wwj_production",       "Exp Centrality (wwj_production); exp(Centrality); Entries", num_bins, *exp_centrality_range)
hist_jet_centrality_wwj_production        = ROOT.TH1F("hist_jet_centrality_wwj_production",        "Jet Centrality (wwj_production); Jet Centrality; Entries", num_bins, *jet_centrality_range)
hist_delta_eta_jj_wwj_production          = ROOT.TH1F("hist_delta_eta_jj_wwj_production",          "Delta Eta jj (wwj_production); Δηjj; Entries", num_bins, *delta_eta_jj_range)
hist_m_w_hadronic_wwj_production          = ROOT.TH1F("hist_m_w_hadronic_wwj_production",          "Hadronic W Mass (wwj_production); m_{W}^{hadronic} [GeV]; Entries", num_bins, *m_w_hadronic_range)
hist_m_w_leptonic_wwj_production          = ROOT.TH1F("hist_m_w_leptonic_wwj_production",          "Leptonic W Mass (wwj_production); m_{W}^{leptonic} [GeV]; Entries", num_bins, *m_w_leptonic_range)
# === Additional Kinematic Observables ===
hist_pt_w_leptonic_wwj_production         = ROOT.TH1F("hist_pt_w_leptonic_wwj_production",         "pT(W^leptonic) (wwj_production); p_{T} [GeV]; Entries", num_bins, 0, 400)
hist_pt_w_hadronic_wwj_production         = ROOT.TH1F("hist_pt_w_hadronic_wwj_production",         "pT(W^hadronic) (wwj_production); p_{T} [GeV]; Entries", num_bins, 0, 400)
hist_delta_phi_lep_met_wwj_production     = ROOT.TH1F("hist_delta_phi_lep_met_wwj_production",     "Δφ(ℓ,MET) (wwj_production); Δφ; Entries", num_bins, 0, np.pi)
hist_mt_w_leptonic_wwj_production         = ROOT.TH1F("hist_mt_w_leptonic_wwj_production",         "MT(W^leptonic) (wwj_production); M_{T} [GeV]; Entries", num_bins, 0, 400)
hist_ht_total_wwj_production              = ROOT.TH1F("hist_ht_total_wwj_production",              "HT Total (wwj_production); H_{T} [GeV]; Entries", num_bins, 0, 1000)
hist_delta_phi_jj_wwj_production          = ROOT.TH1F("hist_delta_phi_jj_wwj_production",          "Δφ(j1, j2) (wwj_production); Δφ; Entries", num_bins, 0, np.pi)
hist_delta_phi_wl_wh_wwj_production       = ROOT.TH1F("hist_delta_phi_wl_wh_wwj_production",       "Δφ(Wlep, Whad) (wwj_production); Δφ; Entries", num_bins, 0, np.pi)
hist_delta_eta_wl_wh_wwj_production       = ROOT.TH1F("hist_delta_eta_wl_wh_wwj_production",       "Δη(Wlep, Whad) (wwj_production); Δη; Entries", num_bins, -6, 6)
hist_m_jj_wwj_production                  = ROOT.TH1F("hist_m_jj_wwj_production",                  "Dijet Invariant Mass (wwj_production); m_{jj} [GeV]; Entries", num_bins, 0, 400)
hist_m_lvjj_wwj_production                = ROOT.TH1F("hist_m_lvjj_wwj_production",                "m(lνjj) (wwj_production); m_{lνjj} [GeV]; Entries", num_bins, 0, 1000)





# === ZZJ PRODUCTION ===
hist_lepton_pt_zzj_production             = ROOT.TH1F("hist_lepton_pt_zzj_production",             "Lepton pT (zzj_production); p_{T} [GeV]; Entries", num_bins, *pt_range_lepton)
hist_leading_jet_pt_zzj_production        = ROOT.TH1F("hist_leading_jet_pt_zzj_production",        "Leading Jet pT (zzj_production); p_{T} [GeV]; Entries", num_bins, *pt_range_jet)
hist_lepton_eta_zzj_production            = ROOT.TH1F("hist_lepton_eta_zzj_production",            "Lepton Eta (zzj_production); #eta; Entries", num_bins, *eta_range)
hist_delta_r_zzj_production               = ROOT.TH1F("hist_delta_r_zzj_production",               "Delta R (zzj_production); ΔR; Entries", num_bins, *delta_r_range)
hist_missing_et_zzj_production            = ROOT.TH1F("hist_missing_et_zzj_production",            "MET (zzj_production); MET [GeV]; Entries", num_bins, *met_range)
hist_subleading_jet_eta_zzj_production    = ROOT.TH1F("hist_subleading_jet_eta_zzj_production",    "Centrality (zzj_production); Centrality; Entries", num_bins, *centrality_range)
hist_leading_jet_eta_zzj_production       = ROOT.TH1F("hist_leading_jet_eta_zzj_production",       "Exp Centrality (zzj_production); exp(Centrality); Entries", num_bins, *exp_centrality_range)
hist_jet_centrality_zzj_production        = ROOT.TH1F("hist_jet_centrality_zzj_production",        "Jet Centrality (zzj_production); Jet Centrality; Entries", num_bins, *jet_centrality_range)
hist_delta_eta_jj_zzj_production          = ROOT.TH1F("hist_delta_eta_jj_zzj_production",          "Delta Eta jj (zzj_production); Δηjj; Entries", num_bins, *delta_eta_jj_range)
hist_m_w_hadronic_zzj_production          = ROOT.TH1F("hist_m_w_hadronic_zzj_production",          "Hadronic W Mass (zzj_production); m_{W}^{hadronic} [GeV]; Entries", num_bins, *m_w_hadronic_range)
hist_m_w_leptonic_zzj_production          = ROOT.TH1F("hist_m_w_leptonic_zzj_production",          "Leptonic W Mass (zzj_production); m_{W}^{leptonic} [GeV]; Entries", num_bins, *m_w_leptonic_range)
# === Additional Observables ===
hist_pt_w_leptonic_zzj_production         = ROOT.TH1F("hist_pt_w_leptonic_zzj_production",         "pT(W^leptonic) (zzj_production); p_{T} [GeV]; Entries", num_bins, 0, 400)
hist_pt_w_hadronic_zzj_production         = ROOT.TH1F("hist_pt_w_hadronic_zzj_production",         "pT(W^hadronic) (zzj_production); p_{T} [GeV]; Entries", num_bins, 0, 400)
hist_delta_phi_lep_met_zzj_production     = ROOT.TH1F("hist_delta_phi_lep_met_zzj_production",     "Δφ(ℓ,MET) (zzj_production); Δφ; Entries", num_bins, 0, np.pi)
hist_mt_w_leptonic_zzj_production         = ROOT.TH1F("hist_mt_w_leptonic_zzj_production",         "MT(W^leptonic) (zzj_production); M_{T} [GeV]; Entries", num_bins, 0, 400)
hist_ht_total_zzj_production              = ROOT.TH1F("hist_ht_total_zzj_production",              "HT Total (zzj_production); H_{T} [GeV]; Entries", num_bins, 0, 1000)
hist_delta_phi_jj_zzj_production          = ROOT.TH1F("hist_delta_phi_jj_zzj_production",          "Δφ(j1, j2) (zzj_production); Δφ; Entries", num_bins, 0, np.pi)
hist_delta_phi_wl_wh_zzj_production       = ROOT.TH1F("hist_delta_phi_wl_wh_zzj_production",       "Δφ(Wlep, Whad) (zzj_production); Δφ; Entries", num_bins, 0, np.pi)
hist_delta_eta_wl_wh_zzj_production       = ROOT.TH1F("hist_delta_eta_wl_wh_zzj_production",       "Δη(Wlep, Whad) (zzj_production); Δη; Entries", num_bins, -6, 6)
hist_m_jj_zzj_production                  = ROOT.TH1F("hist_m_jj_zzj_production",                  "Dijet Invariant Mass (zzj_production); m_{jj} [GeV]; Entries", num_bins, 0, 400)
hist_m_lvjj_zzj_production                = ROOT.TH1F("hist_m_lvjj_zzj_production",                "m(lνjj) (zzj_production); m_{lνjj} [GeV]; Entries", num_bins, 0, 1000)




# === WZJ PRODUCTION ===
hist_lepton_pt_wzj_production             = ROOT.TH1F("hist_lepton_pt_wzj_production",             "Lepton pT (wzj_production); p_{T} [GeV]; Entries", num_bins, *pt_range_lepton)
hist_leading_jet_pt_wzj_production        = ROOT.TH1F("hist_leading_jet_pt_wzj_production",        "Leading Jet pT (wzj_production); p_{T} [GeV]; Entries", num_bins, *pt_range_jet)
hist_lepton_eta_wzj_production            = ROOT.TH1F("hist_lepton_eta_wzj_production",            "Lepton Eta (wzj_production); #eta; Entries", num_bins, *eta_range)
hist_delta_r_wzj_production               = ROOT.TH1F("hist_delta_r_wzj_production",               "Delta R (wzj_production); ΔR; Entries", num_bins, *delta_r_range)
hist_missing_et_wzj_production            = ROOT.TH1F("hist_missing_et_wzj_production",            "MET (wzj_production); MET [GeV]; Entries", num_bins, *met_range)
hist_subleading_jet_eta_wzj_production    = ROOT.TH1F("hist_subleading_jet_eta_wzj_production",    "Centrality (wzj_production); Centrality; Entries", num_bins, *centrality_range)
hist_leading_jet_eta_wzj_production       = ROOT.TH1F("hist_leading_jet_eta_wzj_production",       "Exp Centrality (wzj_production); exp(Centrality); Entries", num_bins, *exp_centrality_range)
hist_jet_centrality_wzj_production        = ROOT.TH1F("hist_jet_centrality_wzj_production",        "Jet Centrality (wzj_production); Jet Centrality; Entries", num_bins, *jet_centrality_range)
hist_delta_eta_jj_wzj_production          = ROOT.TH1F("hist_delta_eta_jj_wzj_production",          "Delta Eta jj (wzj_production); Δηjj; Entries", num_bins, *delta_eta_jj_range)
hist_m_w_hadronic_wzj_production          = ROOT.TH1F("hist_m_w_hadronic_wzj_production",          "Hadronic W Mass (wzj_production); m_{W}^{hadronic} [GeV]; Entries", num_bins, *m_w_hadronic_range)
hist_m_w_leptonic_wzj_production          = ROOT.TH1F("hist_m_w_leptonic_wzj_production",          "Leptonic W Mass (wzj_production); m_{W}^{leptonic} [GeV]; Entries", num_bins, *m_w_leptonic_range)
# === Additional Observables ===
hist_pt_w_leptonic_wzj_production         = ROOT.TH1F("hist_pt_w_leptonic_wzj_production",         "pT(W^leptonic) (wzj_production); p_{T} [GeV]; Entries", num_bins, 0, 400)
hist_pt_w_hadronic_wzj_production         = ROOT.TH1F("hist_pt_w_hadronic_wzj_production",         "pT(W^hadronic) (wzj_production); p_{T} [GeV]; Entries", num_bins, 0, 400)
hist_delta_phi_lep_met_wzj_production     = ROOT.TH1F("hist_delta_phi_lep_met_wzj_production",     "Δφ(ℓ,MET) (wzj_production); Δφ; Entries", num_bins, 0, np.pi)
hist_mt_w_leptonic_wzj_production         = ROOT.TH1F("hist_mt_w_leptonic_wzj_production",         "MT(W^leptonic) (wzj_production); M_{T} [GeV]; Entries", num_bins, 0, 400)
hist_ht_total_wzj_production              = ROOT.TH1F("hist_ht_total_wzj_production",              "HT Total (wzj_production); H_{T} [GeV]; Entries", num_bins, 0, 1000)
hist_delta_phi_jj_wzj_production          = ROOT.TH1F("hist_delta_phi_jj_wzj_production",          "Δφ(j1, j2) (wzj_production); Δφ; Entries", num_bins, 0, np.pi)
hist_delta_phi_wl_wh_wzj_production       = ROOT.TH1F("hist_delta_phi_wl_wh_wzj_production",       "Δφ(Wlep, Whad) (wzj_production); Δφ; Entries", num_bins, 0, np.pi)
hist_delta_eta_wl_wh_wzj_production       = ROOT.TH1F("hist_delta_eta_wl_wh_wzj_production",       "Δη(Wlep, Whad) (wzj_production); Δη; Entries", num_bins, -6, 6)
hist_m_jj_wzj_production                  = ROOT.TH1F("hist_m_jj_wzj_production",                  "Dijet Invariant Mass (wzj_production); m_{jj} [GeV]; Entries", num_bins, 0, 400)
hist_m_lvjj_wzj_production                = ROOT.TH1F("hist_m_lvjj_wzj_production",                "m(lνjj) (wzj_production); m_{lνjj} [GeV]; Entries", num_bins, 0, 1000)







##################################################################


# Dictionary to store efficiencies for each signal
signal_efficiencies = {}

# Process all signal files dynamically
for signal_name, file_path in signal_files.items():
    print(f"Processing signal: {signal_name}")

    # Get the corresponding histograms for this signal
    histograms = signal_histograms[signal_name]

    # Process the signal file and unpack all histograms and efficiencies
    (histograms["hist_lepton_pt"],
     histograms["hist_leading_jet_pt"],
     histograms["hist_lepton_eta"],
     histograms["hist_delta_r"],
     histograms["hist_missing_et"],
     histograms["hist_subleading_jet_eta"],
     histograms["hist_leading_jet_eta"],
     histograms["hist_jet_centrality"],
     histograms["hist_delta_eta_jj"],
     histograms["hist_m_w_leptonic"],
     histograms["hist_m_w_hadronic"],
     histograms["hist_pt_w_leptonic"],
     histograms["hist_pt_w_hadronic"],
     histograms["hist_delta_phi_lep_met"],
     histograms["hist_mt_w_leptonic"],
     histograms["hist_ht_total"],
     histograms["hist_delta_phi_jj"],
     histograms["hist_delta_phi_wl_wh"],
     histograms["hist_delta_eta_wl_wh"],
     histograms["hist_m_jj"],
     histograms["hist_m_lvjj"],
     efficiency_pre_lepton,
     efficiency_pre_jets,
     efficiency_eta_lepton,
     efficiency_jet_centrality,
     efficiency_pre,
     efficiency_final) = process_file(
        file_path,
        histograms["hist_lepton_pt"],
        histograms["hist_leading_jet_pt"],
        histograms["hist_lepton_eta"],
        histograms["hist_delta_r"],
        histograms["hist_missing_et"],
        histograms["hist_subleading_jet_eta"],
        histograms["hist_leading_jet_eta"],
        histograms["hist_jet_centrality"],
        histograms["hist_delta_eta_jj"],
        histograms["hist_m_w_leptonic"],
        histograms["hist_m_w_hadronic"],
        histograms["hist_pt_w_leptonic"],
        histograms["hist_pt_w_hadronic"],
        histograms["hist_delta_phi_lep_met"],
        histograms["hist_mt_w_leptonic"],
        histograms["hist_ht_total"],
        histograms["hist_delta_phi_jj"],
        histograms["hist_delta_phi_wl_wh"],
        histograms["hist_delta_eta_wl_wh"],
        histograms["hist_m_jj"],
        histograms["hist_m_lvjj"]
    )

    # Store efficiencies
    signal_efficiencies[signal_name] = {
        "efficiency_pre_lepton": efficiency_pre_lepton,
        "efficiency_pre_jets": efficiency_pre_jets,
        "efficiency_eta_lepton": efficiency_eta_lepton,
        "efficiency_jet_centrality": efficiency_jet_centrality,
        "efficiency_pre": efficiency_pre,
        "efficiency_final": efficiency_final,
    }



##################################################################


# Auto-generated background process block with full observable set

(hist_lepton_pt_aa_ww, hist_leading_jet_pt_aa_ww, hist_lepton_eta_aa_ww, hist_delta_r_aa_ww,
 hist_missing_et_aa_ww, hist_subleading_jet_eta_aa_ww, hist_leading_jet_eta_aa_ww, hist_jet_centrality_aa_ww,
 hist_delta_eta_jj_aa_ww, hist_m_w_leptonic_aa_ww, hist_m_w_hadronic_aa_ww,
 hist_pt_w_leptonic_aa_ww, hist_pt_w_hadronic_aa_ww, hist_delta_phi_lep_met_aa_ww, hist_mt_w_leptonic_aa_ww,
 hist_ht_total_aa_ww, hist_delta_phi_jj_aa_ww, hist_delta_phi_wl_wh_aa_ww, hist_delta_eta_wl_wh_aa_ww,
 hist_m_jj_aa_ww, hist_m_lvjj_aa_ww,
 background_efficiency_pre_lepton_aa_ww, background_efficiency_pre_jets_aa_ww,
 background_efficiency_eta_lepton_aa_ww, background_efficiency_jet_centrality_aa_ww,
 background_efficiency_pre_aa_ww, background_efficiency_final_aa_ww) = process_file(
    aa_ww_background_file_path,
    hist_lepton_pt_aa_ww,
    hist_leading_jet_pt_aa_ww,
    hist_lepton_eta_aa_ww,
    hist_delta_r_aa_ww,
    hist_missing_et_aa_ww,
    hist_subleading_jet_eta_aa_ww,
    hist_leading_jet_eta_aa_ww,
    hist_jet_centrality_aa_ww,
    hist_delta_eta_jj_aa_ww,
    hist_m_w_leptonic_aa_ww,
    hist_m_w_hadronic_aa_ww,
    hist_pt_w_leptonic_aa_ww,
    hist_pt_w_hadronic_aa_ww,
    hist_delta_phi_lep_met_aa_ww,
    hist_mt_w_leptonic_aa_ww,
    hist_ht_total_aa_ww,
    hist_delta_phi_jj_aa_ww,
    hist_delta_phi_wl_wh_aa_ww,
    hist_delta_eta_wl_wh_aa_ww,
    hist_m_jj_aa_ww,
    hist_m_lvjj_aa_ww
)





(hist_lepton_pt_aa_ttbar, hist_leading_jet_pt_aa_ttbar, hist_lepton_eta_aa_ttbar, hist_delta_r_aa_ttbar,
 hist_missing_et_aa_ttbar, hist_subleading_jet_eta_aa_ttbar, hist_leading_jet_eta_aa_ttbar, hist_jet_centrality_aa_ttbar,
 hist_delta_eta_jj_aa_ttbar, hist_m_w_leptonic_aa_ttbar, hist_m_w_hadronic_aa_ttbar,
 hist_pt_w_leptonic_aa_ttbar, hist_pt_w_hadronic_aa_ttbar, hist_delta_phi_lep_met_aa_ttbar, hist_mt_w_leptonic_aa_ttbar,
 hist_ht_total_aa_ttbar, hist_delta_phi_jj_aa_ttbar, hist_delta_phi_wl_wh_aa_ttbar, hist_delta_eta_wl_wh_aa_ttbar,
 hist_m_jj_aa_ttbar, hist_m_lvjj_aa_ttbar,
 background_efficiency_pre_lepton_aa_ttbar, background_efficiency_pre_jets_aa_ttbar,
 background_efficiency_eta_lepton_aa_ttbar, background_efficiency_jet_centrality_aa_ttbar,
 background_efficiency_pre_aa_ttbar, background_efficiency_final_aa_ttbar) = process_file(
    aa_ttbar_background_file_path,
    hist_lepton_pt_aa_ttbar,
    hist_leading_jet_pt_aa_ttbar,
    hist_lepton_eta_aa_ttbar,
    hist_delta_r_aa_ttbar,
    hist_missing_et_aa_ttbar,
    hist_subleading_jet_eta_aa_ttbar,
    hist_leading_jet_eta_aa_ttbar,
    hist_jet_centrality_aa_ttbar,
    hist_delta_eta_jj_aa_ttbar,
    hist_m_w_leptonic_aa_ttbar,
    hist_m_w_hadronic_aa_ttbar,
    hist_pt_w_leptonic_aa_ttbar,
    hist_pt_w_hadronic_aa_ttbar,
    hist_delta_phi_lep_met_aa_ttbar,
    hist_mt_w_leptonic_aa_ttbar,
    hist_ht_total_aa_ttbar,
    hist_delta_phi_jj_aa_ttbar,
    hist_delta_phi_wl_wh_aa_ttbar,
    hist_delta_eta_wl_wh_aa_ttbar,
    hist_m_jj_aa_ttbar,
    hist_m_lvjj_aa_ttbar
)






(hist_lepton_pt_aa_tautau, hist_leading_jet_pt_aa_tautau, hist_lepton_eta_aa_tautau, hist_delta_r_aa_tautau,
 hist_missing_et_aa_tautau, hist_subleading_jet_eta_aa_tautau, hist_leading_jet_eta_aa_tautau, hist_jet_centrality_aa_tautau,
 hist_delta_eta_jj_aa_tautau, hist_m_w_leptonic_aa_tautau, hist_m_w_hadronic_aa_tautau,
 hist_pt_w_leptonic_aa_tautau, hist_pt_w_hadronic_aa_tautau, hist_delta_phi_lep_met_aa_tautau, hist_mt_w_leptonic_aa_tautau,
 hist_ht_total_aa_tautau, hist_delta_phi_jj_aa_tautau, hist_delta_phi_wl_wh_aa_tautau, hist_delta_eta_wl_wh_aa_tautau,
 hist_m_jj_aa_tautau, hist_m_lvjj_aa_tautau,
 background_efficiency_pre_lepton_aa_tautau, background_efficiency_pre_jets_aa_tautau,
 background_efficiency_eta_lepton_aa_tautau, background_efficiency_jet_centrality_aa_tautau,
 background_efficiency_pre_aa_tautau, background_efficiency_final_aa_tautau) = process_file(
    aa_tautau_background_file_path,
    hist_lepton_pt_aa_tautau,
    hist_leading_jet_pt_aa_tautau,
    hist_lepton_eta_aa_tautau,
    hist_delta_r_aa_tautau,
    hist_missing_et_aa_tautau,
    hist_subleading_jet_eta_aa_tautau,
    hist_leading_jet_eta_aa_tautau,
    hist_jet_centrality_aa_tautau,
    hist_delta_eta_jj_aa_tautau,
    hist_m_w_leptonic_aa_tautau,
    hist_m_w_hadronic_aa_tautau,
    hist_pt_w_leptonic_aa_tautau,
    hist_pt_w_hadronic_aa_tautau,
    hist_delta_phi_lep_met_aa_tautau,
    hist_mt_w_leptonic_aa_tautau,
    hist_ht_total_aa_tautau,
    hist_delta_phi_jj_aa_tautau,
    hist_delta_phi_wl_wh_aa_tautau,
    hist_delta_eta_wl_wh_aa_tautau,
    hist_m_jj_aa_tautau,
    hist_m_lvjj_aa_tautau
)





(hist_lepton_pt_aa_tautau_inel, hist_leading_jet_pt_aa_tautau_inel, hist_lepton_eta_aa_tautau_inel, hist_delta_r_aa_tautau_inel,
 hist_missing_et_aa_tautau_inel, hist_subleading_jet_eta_aa_tautau_inel, hist_leading_jet_eta_aa_tautau_inel, hist_jet_centrality_aa_tautau_inel,
 hist_delta_eta_jj_aa_tautau_inel, hist_m_w_leptonic_aa_tautau_inel, hist_m_w_hadronic_aa_tautau_inel,
 hist_pt_w_leptonic_aa_tautau_inel, hist_pt_w_hadronic_aa_tautau_inel, hist_delta_phi_lep_met_aa_tautau_inel, hist_mt_w_leptonic_aa_tautau_inel,
 hist_ht_total_aa_tautau_inel, hist_delta_phi_jj_aa_tautau_inel, hist_delta_phi_wl_wh_aa_tautau_inel, hist_delta_eta_wl_wh_aa_tautau_inel,
 hist_m_jj_aa_tautau_inel, hist_m_lvjj_aa_tautau_inel,
 background_efficiency_pre_lepton_aa_tautau_inel, background_efficiency_pre_jets_aa_tautau_inel,
 background_efficiency_eta_lepton_aa_tautau_inel, background_efficiency_jet_centrality_aa_tautau_inel,
 background_efficiency_pre_aa_tautau_inel, background_efficiency_final_aa_tautau_inel) = process_file(
    aa_tautau_inel_background_file_path,
    hist_lepton_pt_aa_tautau_inel,
    hist_leading_jet_pt_aa_tautau_inel,
    hist_lepton_eta_aa_tautau_inel,
    hist_delta_r_aa_tautau_inel,
    hist_missing_et_aa_tautau_inel,
    hist_subleading_jet_eta_aa_tautau_inel,
    hist_leading_jet_eta_aa_tautau_inel,
    hist_jet_centrality_aa_tautau_inel,
    hist_delta_eta_jj_aa_tautau_inel,
    hist_m_w_leptonic_aa_tautau_inel,
    hist_m_w_hadronic_aa_tautau_inel,
    hist_pt_w_leptonic_aa_tautau_inel,
    hist_pt_w_hadronic_aa_tautau_inel,
    hist_delta_phi_lep_met_aa_tautau_inel,
    hist_mt_w_leptonic_aa_tautau_inel,
    hist_ht_total_aa_tautau_inel,
    hist_delta_phi_jj_aa_tautau_inel,
    hist_delta_phi_wl_wh_aa_tautau_inel,
    hist_delta_eta_wl_wh_aa_tautau_inel,
    hist_m_jj_aa_tautau_inel,
    hist_m_lvjj_aa_tautau_inel
)






(hist_lepton_pt_inclusive_ttbar, hist_leading_jet_pt_inclusive_ttbar, hist_lepton_eta_inclusive_ttbar, hist_delta_r_inclusive_ttbar,
 hist_missing_et_inclusive_ttbar, hist_subleading_jet_eta_inclusive_ttbar, hist_leading_jet_eta_inclusive_ttbar, hist_jet_centrality_inclusive_ttbar,
 hist_delta_eta_jj_inclusive_ttbar, hist_m_w_leptonic_inclusive_ttbar, hist_m_w_hadronic_inclusive_ttbar,
 hist_pt_w_leptonic_inclusive_ttbar, hist_pt_w_hadronic_inclusive_ttbar, hist_delta_phi_lep_met_inclusive_ttbar, hist_mt_w_leptonic_inclusive_ttbar,
 hist_ht_total_inclusive_ttbar, hist_delta_phi_jj_inclusive_ttbar, hist_delta_phi_wl_wh_inclusive_ttbar, hist_delta_eta_wl_wh_inclusive_ttbar,
 hist_m_jj_inclusive_ttbar, hist_m_lvjj_inclusive_ttbar,
 background_efficiency_pre_lepton_inclusive_ttbar, background_efficiency_pre_jets_inclusive_ttbar,
 background_efficiency_eta_lepton_inclusive_ttbar, background_efficiency_jet_centrality_inclusive_ttbar,
 background_efficiency_pre_inclusive_ttbar, background_efficiency_final_inclusive_ttbar) = process_file(
    inclusive_ttbar_background_file_path,
    hist_lepton_pt_inclusive_ttbar,
    hist_leading_jet_pt_inclusive_ttbar,
    hist_lepton_eta_inclusive_ttbar,
    hist_delta_r_inclusive_ttbar,
    hist_missing_et_inclusive_ttbar,
    hist_subleading_jet_eta_inclusive_ttbar,
    hist_leading_jet_eta_inclusive_ttbar,
    hist_jet_centrality_inclusive_ttbar,
    hist_delta_eta_jj_inclusive_ttbar,
    hist_m_w_leptonic_inclusive_ttbar,
    hist_m_w_hadronic_inclusive_ttbar,
    hist_pt_w_leptonic_inclusive_ttbar,
    hist_pt_w_hadronic_inclusive_ttbar,
    hist_delta_phi_lep_met_inclusive_ttbar,
    hist_mt_w_leptonic_inclusive_ttbar,
    hist_ht_total_inclusive_ttbar,
    hist_delta_phi_jj_inclusive_ttbar,
    hist_delta_phi_wl_wh_inclusive_ttbar,
    hist_delta_eta_wl_wh_inclusive_ttbar,
    hist_m_jj_inclusive_ttbar,
    hist_m_lvjj_inclusive_ttbar
)




(hist_lepton_pt_single_top, hist_leading_jet_pt_single_top, hist_lepton_eta_single_top, hist_delta_r_single_top,
 hist_missing_et_single_top, hist_subleading_jet_eta_single_top, hist_leading_jet_eta_single_top, hist_jet_centrality_single_top,
 hist_delta_eta_jj_single_top, hist_m_w_leptonic_single_top, hist_m_w_hadronic_single_top,
 hist_pt_w_leptonic_single_top, hist_pt_w_hadronic_single_top, hist_delta_phi_lep_met_single_top, hist_mt_w_leptonic_single_top,
 hist_ht_total_single_top, hist_delta_phi_jj_single_top, hist_delta_phi_wl_wh_single_top, hist_delta_eta_wl_wh_single_top,
 hist_m_jj_single_top, hist_m_lvjj_single_top,
 background_efficiency_pre_lepton_single_top, background_efficiency_pre_jets_single_top,
 background_efficiency_eta_lepton_single_top, background_efficiency_jet_centrality_single_top,
 background_efficiency_pre_single_top, background_efficiency_final_single_top) = process_file(
    single_top_background_file_path,
    hist_lepton_pt_single_top,
    hist_leading_jet_pt_single_top,
    hist_lepton_eta_single_top,
    hist_delta_r_single_top,
    hist_missing_et_single_top,
    hist_subleading_jet_eta_single_top,
    hist_leading_jet_eta_single_top,
    hist_jet_centrality_single_top,
    hist_delta_eta_jj_single_top,
    hist_m_w_leptonic_single_top,
    hist_m_w_hadronic_single_top,
    hist_pt_w_leptonic_single_top,
    hist_pt_w_hadronic_single_top,
    hist_delta_phi_lep_met_single_top,
    hist_mt_w_leptonic_single_top,
    hist_ht_total_single_top,
    hist_delta_phi_jj_single_top,
    hist_delta_phi_wl_wh_single_top,
    hist_delta_eta_wl_wh_single_top,
    hist_m_jj_single_top,
    hist_m_lvjj_single_top
)






(hist_lepton_pt_w_production, hist_leading_jet_pt_w_production, hist_lepton_eta_w_production, hist_delta_r_w_production,
 hist_missing_et_w_production, hist_subleading_jet_eta_w_production, hist_leading_jet_eta_w_production, hist_jet_centrality_w_production,
 hist_delta_eta_jj_w_production, hist_m_w_leptonic_w_production, hist_m_w_hadronic_w_production,
 hist_pt_w_leptonic_w_production, hist_pt_w_hadronic_w_production, hist_delta_phi_lep_met_w_production, hist_mt_w_leptonic_w_production,
 hist_ht_total_w_production, hist_delta_phi_jj_w_production, hist_delta_phi_wl_wh_w_production, hist_delta_eta_wl_wh_w_production,
 hist_m_jj_w_production, hist_m_lvjj_w_production,
 background_efficiency_pre_lepton_w_production, background_efficiency_pre_jets_w_production,
 background_efficiency_eta_lepton_w_production, background_efficiency_jet_centrality_w_production,
 background_efficiency_pre_w_production, background_efficiency_final_w_production) = process_file(
    w_production_background_file_path,
    hist_lepton_pt_w_production,
    hist_leading_jet_pt_w_production,
    hist_lepton_eta_w_production,
    hist_delta_r_w_production,
    hist_missing_et_w_production,
    hist_subleading_jet_eta_w_production,
    hist_leading_jet_eta_w_production,
    hist_jet_centrality_w_production,
    hist_delta_eta_jj_w_production,
    hist_m_w_leptonic_w_production,
    hist_m_w_hadronic_w_production,
    hist_pt_w_leptonic_w_production,
    hist_pt_w_hadronic_w_production,
    hist_delta_phi_lep_met_w_production,
    hist_mt_w_leptonic_w_production,
    hist_ht_total_w_production,
    hist_delta_phi_jj_w_production,
    hist_delta_phi_wl_wh_w_production,
    hist_delta_eta_wl_wh_w_production,
    hist_m_jj_w_production,
    hist_m_lvjj_w_production
)






(hist_lepton_pt_z_production, hist_leading_jet_pt_z_production, hist_lepton_eta_z_production, hist_delta_r_z_production,
 hist_missing_et_z_production, hist_subleading_jet_eta_z_production, hist_leading_jet_eta_z_production, hist_jet_centrality_z_production,
 hist_delta_eta_jj_z_production, hist_m_w_leptonic_z_production, hist_m_w_hadronic_z_production,
 hist_pt_w_leptonic_z_production, hist_pt_w_hadronic_z_production, hist_delta_phi_lep_met_z_production, hist_mt_w_leptonic_z_production,
 hist_ht_total_z_production, hist_delta_phi_jj_z_production, hist_delta_phi_wl_wh_z_production, hist_delta_eta_wl_wh_z_production,
 hist_m_jj_z_production, hist_m_lvjj_z_production,
 background_efficiency_pre_lepton_z_production, background_efficiency_pre_jets_z_production,
 background_efficiency_eta_lepton_z_production, background_efficiency_jet_centrality_z_production,
 background_efficiency_pre_z_production, background_efficiency_final_z_production) = process_file(
    z_production_background_file_path,
    hist_lepton_pt_z_production,
    hist_leading_jet_pt_z_production,
    hist_lepton_eta_z_production,
    hist_delta_r_z_production,
    hist_missing_et_z_production,
    hist_subleading_jet_eta_z_production,
    hist_leading_jet_eta_z_production,
    hist_jet_centrality_z_production,
    hist_delta_eta_jj_z_production,
    hist_m_w_leptonic_z_production,
    hist_m_w_hadronic_z_production,
    hist_pt_w_leptonic_z_production,
    hist_pt_w_hadronic_z_production,
    hist_delta_phi_lep_met_z_production,
    hist_mt_w_leptonic_z_production,
    hist_ht_total_z_production,
    hist_delta_phi_jj_z_production,
    hist_delta_phi_wl_wh_z_production,
    hist_delta_eta_wl_wh_z_production,
    hist_m_jj_z_production,
    hist_m_lvjj_z_production
)





(hist_lepton_pt_wwj_production, hist_leading_jet_pt_wwj_production, hist_lepton_eta_wwj_production, hist_delta_r_wwj_production,
 hist_missing_et_wwj_production, hist_subleading_jet_eta_wwj_production, hist_leading_jet_eta_wwj_production, hist_jet_centrality_wwj_production,
 hist_delta_eta_jj_wwj_production, hist_m_w_leptonic_wwj_production, hist_m_w_hadronic_wwj_production,
 hist_pt_w_leptonic_wwj_production, hist_pt_w_hadronic_wwj_production, hist_delta_phi_lep_met_wwj_production, hist_mt_w_leptonic_wwj_production,
 hist_ht_total_wwj_production, hist_delta_phi_jj_wwj_production, hist_delta_phi_wl_wh_wwj_production, hist_delta_eta_wl_wh_wwj_production,
 hist_m_jj_wwj_production, hist_m_lvjj_wwj_production,
 background_efficiency_pre_lepton_wwj_production, background_efficiency_pre_jets_wwj_production,
 background_efficiency_eta_lepton_wwj_production, background_efficiency_jet_centrality_wwj_production,
 background_efficiency_pre_wwj_production, background_efficiency_final_wwj_production) = process_file(
    wwj_production_background_file_path,
    hist_lepton_pt_wwj_production,
    hist_leading_jet_pt_wwj_production,
    hist_lepton_eta_wwj_production,
    hist_delta_r_wwj_production,
    hist_missing_et_wwj_production,
    hist_subleading_jet_eta_wwj_production,
    hist_leading_jet_eta_wwj_production,
    hist_jet_centrality_wwj_production,
    hist_delta_eta_jj_wwj_production,
    hist_m_w_leptonic_wwj_production,
    hist_m_w_hadronic_wwj_production,
    hist_pt_w_leptonic_wwj_production,
    hist_pt_w_hadronic_wwj_production,
    hist_delta_phi_lep_met_wwj_production,
    hist_mt_w_leptonic_wwj_production,
    hist_ht_total_wwj_production,
    hist_delta_phi_jj_wwj_production,
    hist_delta_phi_wl_wh_wwj_production,
    hist_delta_eta_wl_wh_wwj_production,
    hist_m_jj_wwj_production,
    hist_m_lvjj_wwj_production
)




(hist_lepton_pt_zzj_production, hist_leading_jet_pt_zzj_production, hist_lepton_eta_zzj_production, hist_delta_r_zzj_production,
 hist_missing_et_zzj_production, hist_subleading_jet_eta_zzj_production, hist_leading_jet_eta_zzj_production, hist_jet_centrality_zzj_production,
 hist_delta_eta_jj_zzj_production, hist_m_w_leptonic_zzj_production, hist_m_w_hadronic_zzj_production,
 hist_pt_w_leptonic_zzj_production, hist_pt_w_hadronic_zzj_production, hist_delta_phi_lep_met_zzj_production, hist_mt_w_leptonic_zzj_production,
 hist_ht_total_zzj_production, hist_delta_phi_jj_zzj_production, hist_delta_phi_wl_wh_zzj_production, hist_delta_eta_wl_wh_zzj_production,
 hist_m_jj_zzj_production, hist_m_lvjj_zzj_production,
 background_efficiency_pre_lepton_zzj_production, background_efficiency_pre_jets_zzj_production,
 background_efficiency_eta_lepton_zzj_production, background_efficiency_jet_centrality_zzj_production,
 background_efficiency_pre_zzj_production, background_efficiency_final_zzj_production) = process_file(
    zzj_production_background_file_path,
    hist_lepton_pt_zzj_production,
    hist_leading_jet_pt_zzj_production,
    hist_lepton_eta_zzj_production,
    hist_delta_r_zzj_production,
    hist_missing_et_zzj_production,
    hist_subleading_jet_eta_zzj_production,
    hist_leading_jet_eta_zzj_production,
    hist_jet_centrality_zzj_production,
    hist_delta_eta_jj_zzj_production,
    hist_m_w_leptonic_zzj_production,
    hist_m_w_hadronic_zzj_production,
    hist_pt_w_leptonic_zzj_production,
    hist_pt_w_hadronic_zzj_production,
    hist_delta_phi_lep_met_zzj_production,
    hist_mt_w_leptonic_zzj_production,
    hist_ht_total_zzj_production,
    hist_delta_phi_jj_zzj_production,
    hist_delta_phi_wl_wh_zzj_production,
    hist_delta_eta_wl_wh_zzj_production,
    hist_m_jj_zzj_production,
    hist_m_lvjj_zzj_production
)





(hist_lepton_pt_wzj_production, hist_leading_jet_pt_wzj_production, hist_lepton_eta_wzj_production, hist_delta_r_wzj_production,
 hist_missing_et_wzj_production, hist_subleading_jet_eta_wzj_production, hist_leading_jet_eta_wzj_production, hist_jet_centrality_wzj_production,
 hist_delta_eta_jj_wzj_production, hist_m_w_leptonic_wzj_production, hist_m_w_hadronic_wzj_production,
 hist_pt_w_leptonic_wzj_production, hist_pt_w_hadronic_wzj_production, hist_delta_phi_lep_met_wzj_production, hist_mt_w_leptonic_wzj_production,
 hist_ht_total_wzj_production, hist_delta_phi_jj_wzj_production, hist_delta_phi_wl_wh_wzj_production, hist_delta_eta_wl_wh_wzj_production,
 hist_m_jj_wzj_production, hist_m_lvjj_wzj_production,
 background_efficiency_pre_lepton_wzj_production, background_efficiency_pre_jets_wzj_production,
 background_efficiency_eta_lepton_wzj_production, background_efficiency_jet_centrality_wzj_production,
 background_efficiency_pre_wzj_production, background_efficiency_final_wzj_production) = process_file(
    wzj_production_background_file_path,
    hist_lepton_pt_wzj_production,
    hist_leading_jet_pt_wzj_production,
    hist_lepton_eta_wzj_production,
    hist_delta_r_wzj_production,
    hist_missing_et_wzj_production,
    hist_subleading_jet_eta_wzj_production,
    hist_leading_jet_eta_wzj_production,
    hist_jet_centrality_wzj_production,
    hist_delta_eta_jj_wzj_production,
    hist_m_w_leptonic_wzj_production,
    hist_m_w_hadronic_wzj_production,
    hist_pt_w_leptonic_wzj_production,
    hist_pt_w_hadronic_wzj_production,
    hist_delta_phi_lep_met_wzj_production,
    hist_mt_w_leptonic_wzj_production,
    hist_ht_total_wzj_production,
    hist_delta_phi_jj_wzj_production,
    hist_delta_phi_wl_wh_wzj_production,
    hist_delta_eta_wl_wh_wzj_production,
    hist_m_jj_wzj_production,
    hist_m_lvjj_wzj_production
)








####################################################################


print("\n=== Signal Selection Efficiencies ===")
for signal_name, efficiencies in signal_efficiencies.items():
    print(f"{signal_name} Selection Efficiency (1lepton + lepton pT>10 GeV + MET >10 GeV)     : {efficiencies['efficiency_pre_lepton']:.2%}")
    print(f"{signal_name} Selection Efficiency (2jest + iet pT > 10 GeV)        : {efficiencies['efficiency_pre_jets']:.2%}")
    print(f"{signal_name} Selection Efficiency (Lepton η)      : {efficiencies['efficiency_eta_lepton']:.2%}")
    print(f"{signal_name} Selection Efficiency (Jet Centrality) : {efficiencies['efficiency_jet_centrality']:.2%}")
#    print(f"{signal_name} Selection Efficiency (Pre)           : {efficiencies['efficiency_pre']:.2%}")
    print(f"{signal_name} Selection Efficiency (W Mass Window)         : {efficiencies['efficiency_final']:.2%}")
    print("-" * 50)




print("\n=== Background Selection Efficiencies ===")
def print_efficiencies(tag, lep_pt, jet_pt, lep_eta, jet_cent, pre, final):
    print(f"{tag:<15}: 1lepton + lepton pT>10 GeV + MET >10 GeV = {lep_pt:.2%}, 2jest + iet pT > 10 GeV = {jet_pt:.2%}, Lepton η = {lep_eta:.2%}, Jet Centrality = {jet_cent:.2%}, W Mass Window = {final:.2%}")  # Pre = {pre:.2%},

print_efficiencies("aa ww", background_efficiency_pre_lepton_aa_ww, background_efficiency_pre_jets_aa_ww,
                   background_efficiency_eta_lepton_aa_ww, background_efficiency_jet_centrality_aa_ww,
                   background_efficiency_pre_aa_ww, background_efficiency_final_aa_ww)

print_efficiencies("aa ttbar", background_efficiency_pre_lepton_aa_ttbar, background_efficiency_pre_jets_aa_ttbar,
                   background_efficiency_eta_lepton_aa_ttbar, background_efficiency_jet_centrality_aa_ttbar,
                   background_efficiency_pre_aa_ttbar, background_efficiency_final_aa_ttbar)

print_efficiencies("aa tautau", background_efficiency_pre_lepton_aa_tautau, background_efficiency_pre_jets_aa_tautau,
                   background_efficiency_eta_lepton_aa_tautau, background_efficiency_jet_centrality_aa_tautau,
                   background_efficiency_pre_aa_tautau, background_efficiency_final_aa_tautau)

print_efficiencies("aa tautau_inel", background_efficiency_pre_lepton_aa_tautau_inel, background_efficiency_pre_jets_aa_tautau_inel,
                   background_efficiency_eta_lepton_aa_tautau_inel, background_efficiency_jet_centrality_aa_tautau_inel,
                   background_efficiency_pre_aa_tautau_inel, background_efficiency_final_aa_tautau_inel)

print_efficiencies("inclusive ttbar", background_efficiency_pre_lepton_inclusive_ttbar, background_efficiency_pre_jets_inclusive_ttbar,
                   background_efficiency_eta_lepton_inclusive_ttbar, background_efficiency_jet_centrality_inclusive_ttbar,
                   background_efficiency_pre_inclusive_ttbar, background_efficiency_final_inclusive_ttbar)

print_efficiencies("single top", background_efficiency_pre_lepton_single_top, background_efficiency_pre_jets_single_top,
                   background_efficiency_eta_lepton_single_top, background_efficiency_jet_centrality_single_top,
                   background_efficiency_pre_single_top, background_efficiency_final_single_top)

print_efficiencies("w production", background_efficiency_pre_lepton_w_production, background_efficiency_pre_jets_w_production,
                   background_efficiency_eta_lepton_w_production, background_efficiency_jet_centrality_w_production,
                   background_efficiency_pre_w_production, background_efficiency_final_w_production)

print_efficiencies("z production", background_efficiency_pre_lepton_z_production, background_efficiency_pre_jets_z_production,
                   background_efficiency_eta_lepton_z_production, background_efficiency_jet_centrality_z_production,
                   background_efficiency_pre_z_production, background_efficiency_final_z_production)

print_efficiencies("wwj production", background_efficiency_pre_lepton_wwj_production, background_efficiency_pre_jets_wwj_production,
                   background_efficiency_eta_lepton_wwj_production, background_efficiency_jet_centrality_wwj_production,
                   background_efficiency_pre_wwj_production, background_efficiency_final_wwj_production)

print_efficiencies("zzj production", background_efficiency_pre_lepton_zzj_production, background_efficiency_pre_jets_zzj_production,
                   background_efficiency_eta_lepton_zzj_production, background_efficiency_jet_centrality_zzj_production,
                   background_efficiency_pre_zzj_production, background_efficiency_final_zzj_production)

print_efficiencies("wzj production", background_efficiency_pre_lepton_wzj_production, background_efficiency_pre_jets_wzj_production,
                   background_efficiency_eta_lepton_wzj_production, background_efficiency_jet_centrality_wzj_production,
                   background_efficiency_pre_wzj_production, background_efficiency_final_wzj_production)

print("=" * 50)








##############################################################


# Dictionary to store differential cross-sections for each signal
signal_dsigma = {}

# Calculate differential cross-sections for all signals
for signal_name, histograms in signal_histograms.items():
    print(f"Calculating differential cross-sections for {signal_name}...")

    cross_section = signal_cross_sections[signal_name]  # Get cross-section for this signal

    signal_dsigma[signal_name] = {
        "pt_bins_lepton"          : calculate_dsigma(histograms["hist_lepton_pt"], cross_section, bin_width_pt_lepton),
        "pt_bins_jet"             : calculate_dsigma(histograms["hist_leading_jet_pt"], cross_section, bin_width_pt_jet),
        "eta_bins_lepton"         : calculate_dsigma(histograms["hist_lepton_eta"], cross_section, bin_width_eta),
        "delta_r_bins"            : calculate_dsigma(histograms["hist_delta_r"], cross_section, bin_width_delta_r),
        "met_bins"                : calculate_dsigma(histograms["hist_missing_et"], cross_section, bin_width_met),
        "centrality_bins"         : calculate_dsigma(histograms["hist_subleading_jet_eta"], cross_section, bin_width_centrality),
        "exp_centrality_bins"     : calculate_dsigma(histograms["hist_leading_jet_eta"], cross_section, bin_width_exp_centrality),
        "jet_centrality_bins"     : calculate_dsigma(histograms["hist_jet_centrality"], cross_section, bin_width_jet_centrality),
        "delta_eta_jj_bins"       : calculate_dsigma(histograms["hist_delta_eta_jj"], cross_section, bin_width_delta_eta_jj),
        "m_w_hadronic_bins"       : calculate_dsigma(histograms["hist_m_w_hadronic"], cross_section, bin_width_m_w_hadronic),
        "m_w_leptonic_bins"       : calculate_dsigma(histograms["hist_m_w_leptonic"], cross_section, bin_width_m_w_leptonic),

        # ✅ New variables
        "pt_w_leptonic_bins"      : calculate_dsigma(histograms["hist_pt_w_leptonic"], cross_section, bin_width_pt_jet),
        "pt_w_hadronic_bins"      : calculate_dsigma(histograms["hist_pt_w_hadronic"], cross_section, bin_width_pt_jet),
        "delta_phi_lep_met_bins"  : calculate_dsigma(histograms["hist_delta_phi_lep_met"], cross_section, bin_width_delta_r),
        "mt_w_leptonic_bins"      : calculate_dsigma(histograms["hist_mt_w_leptonic"], cross_section, bin_width_m_w_leptonic),
        "ht_total_bins"           : calculate_dsigma(histograms["hist_ht_total"], cross_section, bin_width_pt_jet),
        "delta_phi_jj_bins"       : calculate_dsigma(histograms["hist_delta_phi_jj"], cross_section, bin_width_delta_r),
        "delta_phi_wl_wh_bins"    : calculate_dsigma(histograms["hist_delta_phi_wl_wh"], cross_section, bin_width_delta_r),
        "delta_eta_wl_wh_bins"    : calculate_dsigma(histograms["hist_delta_eta_wl_wh"], cross_section, bin_width_delta_eta_jj),
        "m_jj_bins"               : calculate_dsigma(histograms["hist_m_jj"], cross_section, bin_width_m_w_hadronic),
        "m_lvjj_bins"             : calculate_dsigma(histograms["hist_m_lvjj"], cross_section, bin_width_m_w_leptonic),
    }





# === Differential cross-sections for aa_ww background ===
background_dsigma_aa_ww = {
    "pt_bins_lepton"          : calculate_dsigma(hist_lepton_pt_aa_ww, aa_ww_background_cross_section, bin_width_pt_lepton),
    "pt_bins_jet"             : calculate_dsigma(hist_leading_jet_pt_aa_ww, aa_ww_background_cross_section, bin_width_pt_jet),
    "eta_bins_lepton"         : calculate_dsigma(hist_lepton_eta_aa_ww, aa_ww_background_cross_section, bin_width_eta),
    "delta_r_bins"            : calculate_dsigma(hist_delta_r_aa_ww, aa_ww_background_cross_section, bin_width_delta_r),
    "met_bins"                : calculate_dsigma(hist_missing_et_aa_ww, aa_ww_background_cross_section, bin_width_met),
    "centrality_bins"         : calculate_dsigma(hist_subleading_jet_eta_aa_ww, aa_ww_background_cross_section, bin_width_centrality),
    "exp_centrality_bins"     : calculate_dsigma(hist_leading_jet_eta_aa_ww, aa_ww_background_cross_section, bin_width_exp_centrality),
    "jet_centrality_bins"     : calculate_dsigma(hist_jet_centrality_aa_ww, aa_ww_background_cross_section, bin_width_jet_centrality),
    "delta_eta_jj_bins"       : calculate_dsigma(hist_delta_eta_jj_aa_ww, aa_ww_background_cross_section, bin_width_delta_eta_jj),
    "m_w_hadronic_bins"       : calculate_dsigma(hist_m_w_hadronic_aa_ww, aa_ww_background_cross_section, bin_width_m_w_hadronic),
    "m_w_leptonic_bins"       : calculate_dsigma(hist_m_w_leptonic_aa_ww, aa_ww_background_cross_section, bin_width_m_w_leptonic),

    # ✅ Extended Observables
    "pt_w_leptonic_bins"      : calculate_dsigma(hist_pt_w_leptonic_aa_ww, aa_ww_background_cross_section, bin_width_pt_jet),
    "pt_w_hadronic_bins"      : calculate_dsigma(hist_pt_w_hadronic_aa_ww, aa_ww_background_cross_section, bin_width_pt_jet),
    "delta_phi_lep_met_bins"  : calculate_dsigma(hist_delta_phi_lep_met_aa_ww, aa_ww_background_cross_section, bin_width_delta_r),
    "mt_w_leptonic_bins"      : calculate_dsigma(hist_mt_w_leptonic_aa_ww, aa_ww_background_cross_section, bin_width_m_w_leptonic),
    "ht_total_bins"           : calculate_dsigma(hist_ht_total_aa_ww, aa_ww_background_cross_section, bin_width_pt_jet),
    "delta_phi_jj_bins"       : calculate_dsigma(hist_delta_phi_jj_aa_ww, aa_ww_background_cross_section, bin_width_delta_r),
    "delta_phi_wl_wh_bins"    : calculate_dsigma(hist_delta_phi_wl_wh_aa_ww, aa_ww_background_cross_section, bin_width_delta_r),
    "delta_eta_wl_wh_bins"    : calculate_dsigma(hist_delta_eta_wl_wh_aa_ww, aa_ww_background_cross_section, bin_width_delta_eta_jj),
    "m_jj_bins"               : calculate_dsigma(hist_m_jj_aa_ww, aa_ww_background_cross_section, bin_width_m_w_hadronic),
    "m_lvjj_bins"             : calculate_dsigma(hist_m_lvjj_aa_ww, aa_ww_background_cross_section, bin_width_m_w_leptonic),
}





# === Differential cross-sections for aa_ttbar background ===
background_dsigma_aa_ttbar = {
    "pt_bins_lepton"          : calculate_dsigma(hist_lepton_pt_aa_ttbar, aa_ttbar_background_cross_section, bin_width_pt_lepton),
    "pt_bins_jet"             : calculate_dsigma(hist_leading_jet_pt_aa_ttbar, aa_ttbar_background_cross_section, bin_width_pt_jet),
    "eta_bins_lepton"         : calculate_dsigma(hist_lepton_eta_aa_ttbar, aa_ttbar_background_cross_section, bin_width_eta),
    "delta_r_bins"            : calculate_dsigma(hist_delta_r_aa_ttbar, aa_ttbar_background_cross_section, bin_width_delta_r),
    "met_bins"                : calculate_dsigma(hist_missing_et_aa_ttbar, aa_ttbar_background_cross_section, bin_width_met),
    "centrality_bins"         : calculate_dsigma(hist_subleading_jet_eta_aa_ttbar, aa_ttbar_background_cross_section, bin_width_centrality),
    "exp_centrality_bins"     : calculate_dsigma(hist_leading_jet_eta_aa_ttbar, aa_ttbar_background_cross_section, bin_width_exp_centrality),
    "jet_centrality_bins"     : calculate_dsigma(hist_jet_centrality_aa_ttbar, aa_ttbar_background_cross_section, bin_width_jet_centrality),
    "delta_eta_jj_bins"       : calculate_dsigma(hist_delta_eta_jj_aa_ttbar, aa_ttbar_background_cross_section, bin_width_delta_eta_jj),
    "m_w_hadronic_bins"       : calculate_dsigma(hist_m_w_hadronic_aa_ttbar, aa_ttbar_background_cross_section, bin_width_m_w_hadronic),
    "m_w_leptonic_bins"       : calculate_dsigma(hist_m_w_leptonic_aa_ttbar, aa_ttbar_background_cross_section, bin_width_m_w_leptonic),

    # ✅ Extended Observables
    "pt_w_leptonic_bins"      : calculate_dsigma(hist_pt_w_leptonic_aa_ttbar, aa_ttbar_background_cross_section, bin_width_pt_jet),
    "pt_w_hadronic_bins"      : calculate_dsigma(hist_pt_w_hadronic_aa_ttbar, aa_ttbar_background_cross_section, bin_width_pt_jet),
    "delta_phi_lep_met_bins"  : calculate_dsigma(hist_delta_phi_lep_met_aa_ttbar, aa_ttbar_background_cross_section, bin_width_delta_r),
    "mt_w_leptonic_bins"      : calculate_dsigma(hist_mt_w_leptonic_aa_ttbar, aa_ttbar_background_cross_section, bin_width_m_w_leptonic),
    "ht_total_bins"           : calculate_dsigma(hist_ht_total_aa_ttbar, aa_ttbar_background_cross_section, bin_width_pt_jet),
    "delta_phi_jj_bins"       : calculate_dsigma(hist_delta_phi_jj_aa_ttbar, aa_ttbar_background_cross_section, bin_width_delta_r),
    "delta_phi_wl_wh_bins"    : calculate_dsigma(hist_delta_phi_wl_wh_aa_ttbar, aa_ttbar_background_cross_section, bin_width_delta_r),
    "delta_eta_wl_wh_bins"    : calculate_dsigma(hist_delta_eta_wl_wh_aa_ttbar, aa_ttbar_background_cross_section, bin_width_delta_eta_jj),
    "m_jj_bins"               : calculate_dsigma(hist_m_jj_aa_ttbar, aa_ttbar_background_cross_section, bin_width_m_w_hadronic),
    "m_lvjj_bins"             : calculate_dsigma(hist_m_lvjj_aa_ttbar, aa_ttbar_background_cross_section, bin_width_m_w_leptonic),
}






# === Differential cross-sections for aa_tautau background ===
background_dsigma_aa_tautau = {
    "pt_bins_lepton"          : calculate_dsigma(hist_lepton_pt_aa_tautau, aa_tautau_background_cross_section, bin_width_pt_lepton),
    "pt_bins_jet"             : calculate_dsigma(hist_leading_jet_pt_aa_tautau, aa_tautau_background_cross_section, bin_width_pt_jet),
    "eta_bins_lepton"         : calculate_dsigma(hist_lepton_eta_aa_tautau, aa_tautau_background_cross_section, bin_width_eta),
    "delta_r_bins"            : calculate_dsigma(hist_delta_r_aa_tautau, aa_tautau_background_cross_section, bin_width_delta_r),
    "met_bins"                : calculate_dsigma(hist_missing_et_aa_tautau, aa_tautau_background_cross_section, bin_width_met),
    "centrality_bins"         : calculate_dsigma(hist_subleading_jet_eta_aa_tautau, aa_tautau_background_cross_section, bin_width_centrality),
    "exp_centrality_bins"     : calculate_dsigma(hist_leading_jet_eta_aa_tautau, aa_tautau_background_cross_section, bin_width_exp_centrality),
    "jet_centrality_bins"     : calculate_dsigma(hist_jet_centrality_aa_tautau, aa_tautau_background_cross_section, bin_width_jet_centrality),
    "delta_eta_jj_bins"       : calculate_dsigma(hist_delta_eta_jj_aa_tautau, aa_tautau_background_cross_section, bin_width_delta_eta_jj),
    "m_w_hadronic_bins"       : calculate_dsigma(hist_m_w_hadronic_aa_tautau, aa_tautau_background_cross_section, bin_width_m_w_hadronic),
    "m_w_leptonic_bins"       : calculate_dsigma(hist_m_w_leptonic_aa_tautau, aa_tautau_background_cross_section, bin_width_m_w_leptonic),

    # ✅ Extended observables
    "pt_w_leptonic_bins"      : calculate_dsigma(hist_pt_w_leptonic_aa_tautau, aa_tautau_background_cross_section, bin_width_pt_jet),
    "pt_w_hadronic_bins"      : calculate_dsigma(hist_pt_w_hadronic_aa_tautau, aa_tautau_background_cross_section, bin_width_pt_jet),
    "delta_phi_lep_met_bins"  : calculate_dsigma(hist_delta_phi_lep_met_aa_tautau, aa_tautau_background_cross_section, bin_width_delta_r),
    "mt_w_leptonic_bins"      : calculate_dsigma(hist_mt_w_leptonic_aa_tautau, aa_tautau_background_cross_section, bin_width_m_w_leptonic),
    "ht_total_bins"           : calculate_dsigma(hist_ht_total_aa_tautau, aa_tautau_background_cross_section, bin_width_pt_jet),
    "delta_phi_jj_bins"       : calculate_dsigma(hist_delta_phi_jj_aa_tautau, aa_tautau_background_cross_section, bin_width_delta_r),
    "delta_phi_wl_wh_bins"    : calculate_dsigma(hist_delta_phi_wl_wh_aa_tautau, aa_tautau_background_cross_section, bin_width_delta_r),
    "delta_eta_wl_wh_bins"    : calculate_dsigma(hist_delta_eta_wl_wh_aa_tautau, aa_tautau_background_cross_section, bin_width_delta_eta_jj),
    "m_jj_bins"               : calculate_dsigma(hist_m_jj_aa_tautau, aa_tautau_background_cross_section, bin_width_m_w_hadronic),
    "m_lvjj_bins"             : calculate_dsigma(hist_m_lvjj_aa_tautau, aa_tautau_background_cross_section, bin_width_m_w_leptonic),
}





# === Differential cross-sections for aa_tautau_inel background ===
background_dsigma_aa_tautau_inel = {
    "pt_bins_lepton"          : calculate_dsigma(hist_lepton_pt_aa_tautau_inel, aa_tautau_inel_background_cross_section, bin_width_pt_lepton),
    "pt_bins_jet"             : calculate_dsigma(hist_leading_jet_pt_aa_tautau_inel, aa_tautau_inel_background_cross_section, bin_width_pt_jet),
    "eta_bins_lepton"         : calculate_dsigma(hist_lepton_eta_aa_tautau_inel, aa_tautau_inel_background_cross_section, bin_width_eta),
    "delta_r_bins"            : calculate_dsigma(hist_delta_r_aa_tautau_inel, aa_tautau_inel_background_cross_section, bin_width_delta_r),
    "met_bins"                : calculate_dsigma(hist_missing_et_aa_tautau_inel, aa_tautau_inel_background_cross_section, bin_width_met),
    "centrality_bins"         : calculate_dsigma(hist_subleading_jet_eta_aa_tautau_inel, aa_tautau_inel_background_cross_section, bin_width_centrality),
    "exp_centrality_bins"     : calculate_dsigma(hist_leading_jet_eta_aa_tautau_inel, aa_tautau_inel_background_cross_section, bin_width_exp_centrality),
    "jet_centrality_bins"     : calculate_dsigma(hist_jet_centrality_aa_tautau_inel, aa_tautau_inel_background_cross_section, bin_width_jet_centrality),
    "delta_eta_jj_bins"       : calculate_dsigma(hist_delta_eta_jj_aa_tautau_inel, aa_tautau_inel_background_cross_section, bin_width_delta_eta_jj),
    "m_w_hadronic_bins"       : calculate_dsigma(hist_m_w_hadronic_aa_tautau_inel, aa_tautau_inel_background_cross_section, bin_width_m_w_hadronic),
    "m_w_leptonic_bins"       : calculate_dsigma(hist_m_w_leptonic_aa_tautau_inel, aa_tautau_inel_background_cross_section, bin_width_m_w_leptonic),

    # ✅ Extended observables
    "pt_w_leptonic_bins"      : calculate_dsigma(hist_pt_w_leptonic_aa_tautau_inel, aa_tautau_inel_background_cross_section, bin_width_pt_jet),
    "pt_w_hadronic_bins"      : calculate_dsigma(hist_pt_w_hadronic_aa_tautau_inel, aa_tautau_inel_background_cross_section, bin_width_pt_jet),
    "delta_phi_lep_met_bins"  : calculate_dsigma(hist_delta_phi_lep_met_aa_tautau_inel, aa_tautau_inel_background_cross_section, bin_width_delta_r),
    "mt_w_leptonic_bins"      : calculate_dsigma(hist_mt_w_leptonic_aa_tautau_inel, aa_tautau_inel_background_cross_section, bin_width_m_w_leptonic),
    "ht_total_bins"           : calculate_dsigma(hist_ht_total_aa_tautau_inel, aa_tautau_inel_background_cross_section, bin_width_pt_jet),
    "delta_phi_jj_bins"       : calculate_dsigma(hist_delta_phi_jj_aa_tautau_inel, aa_tautau_inel_background_cross_section, bin_width_delta_r),
    "delta_phi_wl_wh_bins"    : calculate_dsigma(hist_delta_phi_wl_wh_aa_tautau_inel, aa_tautau_inel_background_cross_section, bin_width_delta_r),
    "delta_eta_wl_wh_bins"    : calculate_dsigma(hist_delta_eta_wl_wh_aa_tautau_inel, aa_tautau_inel_background_cross_section, bin_width_delta_eta_jj),
    "m_jj_bins"               : calculate_dsigma(hist_m_jj_aa_tautau_inel, aa_tautau_inel_background_cross_section, bin_width_m_w_hadronic),
    "m_lvjj_bins"             : calculate_dsigma(hist_m_lvjj_aa_tautau_inel, aa_tautau_inel_background_cross_section, bin_width_m_w_leptonic),
}






# === Differential cross-sections for inclusive_ttbar background ===
background_dsigma_inclusive_ttbar = {
    "pt_bins_lepton"          : calculate_dsigma(hist_lepton_pt_inclusive_ttbar, inclusive_ttbar_background_cross_section, bin_width_pt_lepton),
    "pt_bins_jet"             : calculate_dsigma(hist_leading_jet_pt_inclusive_ttbar, inclusive_ttbar_background_cross_section, bin_width_pt_jet),
    "eta_bins_lepton"         : calculate_dsigma(hist_lepton_eta_inclusive_ttbar, inclusive_ttbar_background_cross_section, bin_width_eta),
    "delta_r_bins"            : calculate_dsigma(hist_delta_r_inclusive_ttbar, inclusive_ttbar_background_cross_section, bin_width_delta_r),
    "met_bins"                : calculate_dsigma(hist_missing_et_inclusive_ttbar, inclusive_ttbar_background_cross_section, bin_width_met),
    "centrality_bins"         : calculate_dsigma(hist_subleading_jet_eta_inclusive_ttbar, inclusive_ttbar_background_cross_section, bin_width_centrality),
    "exp_centrality_bins"     : calculate_dsigma(hist_leading_jet_eta_inclusive_ttbar, inclusive_ttbar_background_cross_section, bin_width_exp_centrality),
    "jet_centrality_bins"     : calculate_dsigma(hist_jet_centrality_inclusive_ttbar, inclusive_ttbar_background_cross_section, bin_width_jet_centrality),
    "delta_eta_jj_bins"       : calculate_dsigma(hist_delta_eta_jj_inclusive_ttbar, inclusive_ttbar_background_cross_section, bin_width_delta_eta_jj),
    "m_w_hadronic_bins"       : calculate_dsigma(hist_m_w_hadronic_inclusive_ttbar, inclusive_ttbar_background_cross_section, bin_width_m_w_hadronic),
    "m_w_leptonic_bins"       : calculate_dsigma(hist_m_w_leptonic_inclusive_ttbar, inclusive_ttbar_background_cross_section, bin_width_m_w_leptonic),

    # ✅ Extended observables
    "pt_w_leptonic_bins"      : calculate_dsigma(hist_pt_w_leptonic_inclusive_ttbar, inclusive_ttbar_background_cross_section, bin_width_pt_jet),
    "pt_w_hadronic_bins"      : calculate_dsigma(hist_pt_w_hadronic_inclusive_ttbar, inclusive_ttbar_background_cross_section, bin_width_pt_jet),
    "delta_phi_lep_met_bins"  : calculate_dsigma(hist_delta_phi_lep_met_inclusive_ttbar, inclusive_ttbar_background_cross_section, bin_width_delta_r),
    "mt_w_leptonic_bins"      : calculate_dsigma(hist_mt_w_leptonic_inclusive_ttbar, inclusive_ttbar_background_cross_section, bin_width_m_w_leptonic),
    "ht_total_bins"           : calculate_dsigma(hist_ht_total_inclusive_ttbar, inclusive_ttbar_background_cross_section, bin_width_pt_jet),
    "delta_phi_jj_bins"       : calculate_dsigma(hist_delta_phi_jj_inclusive_ttbar, inclusive_ttbar_background_cross_section, bin_width_delta_r),
    "delta_phi_wl_wh_bins"    : calculate_dsigma(hist_delta_phi_wl_wh_inclusive_ttbar, inclusive_ttbar_background_cross_section, bin_width_delta_r),
    "delta_eta_wl_wh_bins"    : calculate_dsigma(hist_delta_eta_wl_wh_inclusive_ttbar, inclusive_ttbar_background_cross_section, bin_width_delta_eta_jj),
    "m_jj_bins"               : calculate_dsigma(hist_m_jj_inclusive_ttbar, inclusive_ttbar_background_cross_section, bin_width_m_w_hadronic),
    "m_lvjj_bins"             : calculate_dsigma(hist_m_lvjj_inclusive_ttbar, inclusive_ttbar_background_cross_section, bin_width_m_w_leptonic),
}






# === Differential cross-sections for single_top background ===
background_dsigma_single_top = {
    "pt_bins_lepton"         : calculate_dsigma(hist_lepton_pt_single_top, single_top_background_cross_section, bin_width_pt_lepton),
    "pt_bins_jet"            : calculate_dsigma(hist_leading_jet_pt_single_top, single_top_background_cross_section, bin_width_pt_jet),
    "eta_bins_lepton"        : calculate_dsigma(hist_lepton_eta_single_top, single_top_background_cross_section, bin_width_eta),
    "delta_r_bins"           : calculate_dsigma(hist_delta_r_single_top, single_top_background_cross_section, bin_width_delta_r),
    "met_bins"               : calculate_dsigma(hist_missing_et_single_top, single_top_background_cross_section, bin_width_met),
    "centrality_bins"        : calculate_dsigma(hist_subleading_jet_eta_single_top, single_top_background_cross_section, bin_width_centrality),
    "exp_centrality_bins"    : calculate_dsigma(hist_leading_jet_eta_single_top, single_top_background_cross_section, bin_width_exp_centrality),
    "jet_centrality_bins"    : calculate_dsigma(hist_jet_centrality_single_top, single_top_background_cross_section, bin_width_jet_centrality),
    "delta_eta_jj_bins"      : calculate_dsigma(hist_delta_eta_jj_single_top, single_top_background_cross_section, bin_width_delta_eta_jj),
    "m_w_hadronic_bins"      : calculate_dsigma(hist_m_w_hadronic_single_top, single_top_background_cross_section, bin_width_m_w_hadronic),
    "m_w_leptonic_bins"      : calculate_dsigma(hist_m_w_leptonic_single_top, single_top_background_cross_section, bin_width_m_w_leptonic),

    # ✅ Extended observables
    "pt_w_leptonic_bins"     : calculate_dsigma(hist_pt_w_leptonic_single_top, single_top_background_cross_section, bin_width_pt_jet),
    "pt_w_hadronic_bins"     : calculate_dsigma(hist_pt_w_hadronic_single_top, single_top_background_cross_section, bin_width_pt_jet),
    "delta_phi_lep_met_bins" : calculate_dsigma(hist_delta_phi_lep_met_single_top, single_top_background_cross_section, bin_width_delta_r),
    "mt_w_leptonic_bins"     : calculate_dsigma(hist_mt_w_leptonic_single_top, single_top_background_cross_section, bin_width_m_w_leptonic),
    "ht_total_bins"          : calculate_dsigma(hist_ht_total_single_top, single_top_background_cross_section, bin_width_pt_jet),
    "delta_phi_jj_bins"      : calculate_dsigma(hist_delta_phi_jj_single_top, single_top_background_cross_section, bin_width_delta_r),
    "delta_phi_wl_wh_bins"   : calculate_dsigma(hist_delta_phi_wl_wh_single_top, single_top_background_cross_section, bin_width_delta_r),
    "delta_eta_wl_wh_bins"   : calculate_dsigma(hist_delta_eta_wl_wh_single_top, single_top_background_cross_section, bin_width_delta_eta_jj),
    "m_jj_bins"              : calculate_dsigma(hist_m_jj_single_top, single_top_background_cross_section, bin_width_m_w_hadronic),
    "m_lvjj_bins"            : calculate_dsigma(hist_m_lvjj_single_top, single_top_background_cross_section, bin_width_m_w_leptonic),
}






# === Differential cross-sections for w_production background ===
background_dsigma_w_production = {
    "pt_bins_lepton"         : calculate_dsigma(hist_lepton_pt_w_production, w_production_background_cross_section, bin_width_pt_lepton),
    "pt_bins_jet"            : calculate_dsigma(hist_leading_jet_pt_w_production, w_production_background_cross_section, bin_width_pt_jet),
    "eta_bins_lepton"        : calculate_dsigma(hist_lepton_eta_w_production, w_production_background_cross_section, bin_width_eta),
    "delta_r_bins"           : calculate_dsigma(hist_delta_r_w_production, w_production_background_cross_section, bin_width_delta_r),
    "met_bins"               : calculate_dsigma(hist_missing_et_w_production, w_production_background_cross_section, bin_width_met),
    "centrality_bins"        : calculate_dsigma(hist_subleading_jet_eta_w_production, w_production_background_cross_section, bin_width_centrality),
    "exp_centrality_bins"    : calculate_dsigma(hist_leading_jet_eta_w_production, w_production_background_cross_section, bin_width_exp_centrality),
    "jet_centrality_bins"    : calculate_dsigma(hist_jet_centrality_w_production, w_production_background_cross_section, bin_width_jet_centrality),
    "delta_eta_jj_bins"      : calculate_dsigma(hist_delta_eta_jj_w_production, w_production_background_cross_section, bin_width_delta_eta_jj),
    "m_w_hadronic_bins"      : calculate_dsigma(hist_m_w_hadronic_w_production, w_production_background_cross_section, bin_width_m_w_hadronic),
    "m_w_leptonic_bins"      : calculate_dsigma(hist_m_w_leptonic_w_production, w_production_background_cross_section, bin_width_m_w_leptonic),

    # ✅ Extended observables
    "pt_w_leptonic_bins"     : calculate_dsigma(hist_pt_w_leptonic_w_production, w_production_background_cross_section, bin_width_pt_jet),
    "pt_w_hadronic_bins"     : calculate_dsigma(hist_pt_w_hadronic_w_production, w_production_background_cross_section, bin_width_pt_jet),
    "delta_phi_lep_met_bins" : calculate_dsigma(hist_delta_phi_lep_met_w_production, w_production_background_cross_section, bin_width_delta_r),
    "mt_w_leptonic_bins"     : calculate_dsigma(hist_mt_w_leptonic_w_production, w_production_background_cross_section, bin_width_m_w_leptonic),
    "ht_total_bins"          : calculate_dsigma(hist_ht_total_w_production, w_production_background_cross_section, bin_width_pt_jet),
    "delta_phi_jj_bins"      : calculate_dsigma(hist_delta_phi_jj_w_production, w_production_background_cross_section, bin_width_delta_r),
    "delta_phi_wl_wh_bins"   : calculate_dsigma(hist_delta_phi_wl_wh_w_production, w_production_background_cross_section, bin_width_delta_r),
    "delta_eta_wl_wh_bins"   : calculate_dsigma(hist_delta_eta_wl_wh_w_production, w_production_background_cross_section, bin_width_delta_eta_jj),
    "m_jj_bins"              : calculate_dsigma(hist_m_jj_w_production, w_production_background_cross_section, bin_width_m_w_hadronic),
    "m_lvjj_bins"            : calculate_dsigma(hist_m_lvjj_w_production, w_production_background_cross_section, bin_width_m_w_leptonic),
}






# === Differential cross-sections for z_production background ===
background_dsigma_z_production = {
    "pt_bins_lepton"         : calculate_dsigma(hist_lepton_pt_z_production, z_production_background_cross_section, bin_width_pt_lepton),
    "pt_bins_jet"            : calculate_dsigma(hist_leading_jet_pt_z_production, z_production_background_cross_section, bin_width_pt_jet),
    "eta_bins_lepton"        : calculate_dsigma(hist_lepton_eta_z_production, z_production_background_cross_section, bin_width_eta),
    "delta_r_bins"           : calculate_dsigma(hist_delta_r_z_production, z_production_background_cross_section, bin_width_delta_r),
    "met_bins"               : calculate_dsigma(hist_missing_et_z_production, z_production_background_cross_section, bin_width_met),
    "centrality_bins"        : calculate_dsigma(hist_subleading_jet_eta_z_production, z_production_background_cross_section, bin_width_centrality),
    "exp_centrality_bins"    : calculate_dsigma(hist_leading_jet_eta_z_production, z_production_background_cross_section, bin_width_exp_centrality),
    "jet_centrality_bins"    : calculate_dsigma(hist_jet_centrality_z_production, z_production_background_cross_section, bin_width_jet_centrality),
    "delta_eta_jj_bins"      : calculate_dsigma(hist_delta_eta_jj_z_production, z_production_background_cross_section, bin_width_delta_eta_jj),
    "m_w_hadronic_bins"      : calculate_dsigma(hist_m_w_hadronic_z_production, z_production_background_cross_section, bin_width_m_w_hadronic),
    "m_w_leptonic_bins"      : calculate_dsigma(hist_m_w_leptonic_z_production, z_production_background_cross_section, bin_width_m_w_leptonic),

    # ✅ Extended observables
    "pt_w_leptonic_bins"     : calculate_dsigma(hist_pt_w_leptonic_z_production, z_production_background_cross_section, bin_width_pt_jet),
    "pt_w_hadronic_bins"     : calculate_dsigma(hist_pt_w_hadronic_z_production, z_production_background_cross_section, bin_width_pt_jet),
    "delta_phi_lep_met_bins" : calculate_dsigma(hist_delta_phi_lep_met_z_production, z_production_background_cross_section, bin_width_delta_r),
    "mt_w_leptonic_bins"     : calculate_dsigma(hist_mt_w_leptonic_z_production, z_production_background_cross_section, bin_width_m_w_leptonic),
    "ht_total_bins"          : calculate_dsigma(hist_ht_total_z_production, z_production_background_cross_section, bin_width_pt_jet),
    "delta_phi_jj_bins"      : calculate_dsigma(hist_delta_phi_jj_z_production, z_production_background_cross_section, bin_width_delta_r),
    "delta_phi_wl_wh_bins"   : calculate_dsigma(hist_delta_phi_wl_wh_z_production, z_production_background_cross_section, bin_width_delta_r),
    "delta_eta_wl_wh_bins"   : calculate_dsigma(hist_delta_eta_wl_wh_z_production, z_production_background_cross_section, bin_width_delta_eta_jj),
    "m_jj_bins"              : calculate_dsigma(hist_m_jj_z_production, z_production_background_cross_section, bin_width_m_w_hadronic),
    "m_lvjj_bins"            : calculate_dsigma(hist_m_lvjj_z_production, z_production_background_cross_section, bin_width_m_w_leptonic),
}







# === Differential cross-sections for wwj_production background ===
background_dsigma_wwj_production = {
    "pt_bins_lepton"         : calculate_dsigma(hist_lepton_pt_wwj_production, wwj_production_background_cross_section, bin_width_pt_lepton),
    "pt_bins_jet"            : calculate_dsigma(hist_leading_jet_pt_wwj_production, wwj_production_background_cross_section, bin_width_pt_jet),
    "eta_bins_lepton"        : calculate_dsigma(hist_lepton_eta_wwj_production, wwj_production_background_cross_section, bin_width_eta),
    "delta_r_bins"           : calculate_dsigma(hist_delta_r_wwj_production, wwj_production_background_cross_section, bin_width_delta_r),
    "met_bins"               : calculate_dsigma(hist_missing_et_wwj_production, wwj_production_background_cross_section, bin_width_met),
    "centrality_bins"        : calculate_dsigma(hist_subleading_jet_eta_wwj_production, wwj_production_background_cross_section, bin_width_centrality),
    "exp_centrality_bins"    : calculate_dsigma(hist_leading_jet_eta_wwj_production, wwj_production_background_cross_section, bin_width_exp_centrality),
    "jet_centrality_bins"    : calculate_dsigma(hist_jet_centrality_wwj_production, wwj_production_background_cross_section, bin_width_jet_centrality),
    "delta_eta_jj_bins"      : calculate_dsigma(hist_delta_eta_jj_wwj_production, wwj_production_background_cross_section, bin_width_delta_eta_jj),
    "m_w_hadronic_bins"      : calculate_dsigma(hist_m_w_hadronic_wwj_production, wwj_production_background_cross_section, bin_width_m_w_hadronic),
    "m_w_leptonic_bins"      : calculate_dsigma(hist_m_w_leptonic_wwj_production, wwj_production_background_cross_section, bin_width_m_w_leptonic),

    # ✅ Extended observables
    "pt_w_leptonic_bins"     : calculate_dsigma(hist_pt_w_leptonic_wwj_production, wwj_production_background_cross_section, bin_width_pt_jet),
    "pt_w_hadronic_bins"     : calculate_dsigma(hist_pt_w_hadronic_wwj_production, wwj_production_background_cross_section, bin_width_pt_jet),
    "delta_phi_lep_met_bins" : calculate_dsigma(hist_delta_phi_lep_met_wwj_production, wwj_production_background_cross_section, bin_width_delta_r),
    "mt_w_leptonic_bins"     : calculate_dsigma(hist_mt_w_leptonic_wwj_production, wwj_production_background_cross_section, bin_width_m_w_leptonic),
    "ht_total_bins"          : calculate_dsigma(hist_ht_total_wwj_production, wwj_production_background_cross_section, bin_width_pt_jet),
    "delta_phi_jj_bins"      : calculate_dsigma(hist_delta_phi_jj_wwj_production, wwj_production_background_cross_section, bin_width_delta_r),
    "delta_phi_wl_wh_bins"   : calculate_dsigma(hist_delta_phi_wl_wh_wwj_production, wwj_production_background_cross_section, bin_width_delta_r),
    "delta_eta_wl_wh_bins"   : calculate_dsigma(hist_delta_eta_wl_wh_wwj_production, wwj_production_background_cross_section, bin_width_delta_eta_jj),
    "m_jj_bins"              : calculate_dsigma(hist_m_jj_wwj_production, wwj_production_background_cross_section, bin_width_m_w_hadronic),
    "m_lvjj_bins"            : calculate_dsigma(hist_m_lvjj_wwj_production, wwj_production_background_cross_section, bin_width_m_w_leptonic),
}







# === Differential cross-sections for zzj_production background ===
background_dsigma_zzj_production = {
    "pt_bins_lepton"         : calculate_dsigma(hist_lepton_pt_zzj_production, zzj_production_background_cross_section, bin_width_pt_lepton),
    "pt_bins_jet"            : calculate_dsigma(hist_leading_jet_pt_zzj_production, zzj_production_background_cross_section, bin_width_pt_jet),
    "eta_bins_lepton"        : calculate_dsigma(hist_lepton_eta_zzj_production, zzj_production_background_cross_section, bin_width_eta),
    "delta_r_bins"           : calculate_dsigma(hist_delta_r_zzj_production, zzj_production_background_cross_section, bin_width_delta_r),
    "met_bins"               : calculate_dsigma(hist_missing_et_zzj_production, zzj_production_background_cross_section, bin_width_met),
    "centrality_bins"        : calculate_dsigma(hist_subleading_jet_eta_zzj_production, zzj_production_background_cross_section, bin_width_centrality),
    "exp_centrality_bins"    : calculate_dsigma(hist_leading_jet_eta_zzj_production, zzj_production_background_cross_section, bin_width_exp_centrality),
    "jet_centrality_bins"    : calculate_dsigma(hist_jet_centrality_zzj_production, zzj_production_background_cross_section, bin_width_jet_centrality),
    "delta_eta_jj_bins"      : calculate_dsigma(hist_delta_eta_jj_zzj_production, zzj_production_background_cross_section, bin_width_delta_eta_jj),
    "m_w_hadronic_bins"      : calculate_dsigma(hist_m_w_hadronic_zzj_production, zzj_production_background_cross_section, bin_width_m_w_hadronic),
    "m_w_leptonic_bins"      : calculate_dsigma(hist_m_w_leptonic_zzj_production, zzj_production_background_cross_section, bin_width_m_w_leptonic),

    # ✅ Extended observables
    "pt_w_leptonic_bins"     : calculate_dsigma(hist_pt_w_leptonic_zzj_production, zzj_production_background_cross_section, bin_width_pt_jet),
    "pt_w_hadronic_bins"     : calculate_dsigma(hist_pt_w_hadronic_zzj_production, zzj_production_background_cross_section, bin_width_pt_jet),
    "delta_phi_lep_met_bins" : calculate_dsigma(hist_delta_phi_lep_met_zzj_production, zzj_production_background_cross_section, bin_width_delta_r),
    "mt_w_leptonic_bins"     : calculate_dsigma(hist_mt_w_leptonic_zzj_production, zzj_production_background_cross_section, bin_width_m_w_leptonic),
    "ht_total_bins"          : calculate_dsigma(hist_ht_total_zzj_production, zzj_production_background_cross_section, bin_width_pt_jet),
    "delta_phi_jj_bins"      : calculate_dsigma(hist_delta_phi_jj_zzj_production, zzj_production_background_cross_section, bin_width_delta_r),
    "delta_phi_wl_wh_bins"   : calculate_dsigma(hist_delta_phi_wl_wh_zzj_production, zzj_production_background_cross_section, bin_width_delta_r),
    "delta_eta_wl_wh_bins"   : calculate_dsigma(hist_delta_eta_wl_wh_zzj_production, zzj_production_background_cross_section, bin_width_delta_eta_jj),
    "m_jj_bins"              : calculate_dsigma(hist_m_jj_zzj_production, zzj_production_background_cross_section, bin_width_m_w_hadronic),
    "m_lvjj_bins"            : calculate_dsigma(hist_m_lvjj_zzj_production, zzj_production_background_cross_section, bin_width_m_w_leptonic),
}





# === Differential cross-sections for wzj_production background ===
background_dsigma_wzj_production = {
    "pt_bins_lepton"         : calculate_dsigma(hist_lepton_pt_wzj_production, wzj_production_background_cross_section, bin_width_pt_lepton),
    "pt_bins_jet"            : calculate_dsigma(hist_leading_jet_pt_wzj_production, wzj_production_background_cross_section, bin_width_pt_jet),
    "eta_bins_lepton"        : calculate_dsigma(hist_lepton_eta_wzj_production, wzj_production_background_cross_section, bin_width_eta),
    "delta_r_bins"           : calculate_dsigma(hist_delta_r_wzj_production, wzj_production_background_cross_section, bin_width_delta_r),
    "met_bins"               : calculate_dsigma(hist_missing_et_wzj_production, wzj_production_background_cross_section, bin_width_met),
    "centrality_bins"        : calculate_dsigma(hist_subleading_jet_eta_wzj_production, wzj_production_background_cross_section, bin_width_centrality),
    "exp_centrality_bins"    : calculate_dsigma(hist_leading_jet_eta_wzj_production, wzj_production_background_cross_section, bin_width_exp_centrality),
    "jet_centrality_bins"    : calculate_dsigma(hist_jet_centrality_wzj_production, wzj_production_background_cross_section, bin_width_jet_centrality),
    "delta_eta_jj_bins"      : calculate_dsigma(hist_delta_eta_jj_wzj_production, wzj_production_background_cross_section, bin_width_delta_eta_jj),
    "m_w_hadronic_bins"      : calculate_dsigma(hist_m_w_hadronic_wzj_production, wzj_production_background_cross_section, bin_width_m_w_hadronic),
    "m_w_leptonic_bins"      : calculate_dsigma(hist_m_w_leptonic_wzj_production, wzj_production_background_cross_section, bin_width_m_w_leptonic),

    # ✅ Extended observables
    "pt_w_leptonic_bins"     : calculate_dsigma(hist_pt_w_leptonic_wzj_production, wzj_production_background_cross_section, bin_width_pt_jet),
    "pt_w_hadronic_bins"     : calculate_dsigma(hist_pt_w_hadronic_wzj_production, wzj_production_background_cross_section, bin_width_pt_jet),
    "delta_phi_lep_met_bins" : calculate_dsigma(hist_delta_phi_lep_met_wzj_production, wzj_production_background_cross_section, bin_width_delta_r),
    "mt_w_leptonic_bins"     : calculate_dsigma(hist_mt_w_leptonic_wzj_production, wzj_production_background_cross_section, bin_width_m_w_leptonic),
    "ht_total_bins"          : calculate_dsigma(hist_ht_total_wzj_production, wzj_production_background_cross_section, bin_width_pt_jet),
    "delta_phi_jj_bins"      : calculate_dsigma(hist_delta_phi_jj_wzj_production, wzj_production_background_cross_section, bin_width_delta_r),
    "delta_phi_wl_wh_bins"   : calculate_dsigma(hist_delta_phi_wl_wh_wzj_production, wzj_production_background_cross_section, bin_width_delta_r),
    "delta_eta_wl_wh_bins"   : calculate_dsigma(hist_delta_eta_wl_wh_wzj_production, wzj_production_background_cross_section, bin_width_delta_eta_jj),
    "m_jj_bins"              : calculate_dsigma(hist_m_jj_wzj_production, wzj_production_background_cross_section, bin_width_m_w_hadronic),
    "m_lvjj_bins"            : calculate_dsigma(hist_m_lvjj_wzj_production, wzj_production_background_cross_section, bin_width_m_w_leptonic),
}








#=========================================================================
#=========================================================================




# Define colors for each signal
signal_colors = {
    "$FM_{0} / \Lambda^4$": "green",
    "$FM_{1} / \Lambda^4$": "purple",
    "FM2_Lambda4": "red",
    "$FM_{3} / \Lambda^4$": "orange"
}





# ===================================================
# ✅ Lepton \( p_T \) Differential Cross-Section (Semi-Leptonic)
# ===================================================


plt.figure(figsize=(11, 12))
plt.subplots_adjust(left=0.15, right=0.95, bottom=0.12, top=0.95)

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
    (r"$\gamma\gamma \to \tau^+\tau^- (inel)$", background_dsigma_aa_tautau_inel, "brown", (0, (3, 1, 1, 1))),
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
plt.ylim(1e-5, 1e0)
# ✅ Legend, Grid, and Save
plt.legend()
plt.grid(True, linestyle="--", alpha=0.6)
plt.tight_layout()
plt.savefig("/home/hamzeh-khanpour/Documents/GitHub/LHeC_Fast_Simulation/Pythia8_Delphes_semi_leptonic_allsignal_bkgs/differential_cross_section_lepton_pt_allFMsignal_allbkgs_all_decay_mode.pdf", dpi=600)

plt.show()








# ===================================================
# ✅ Lepton \( p_T \) Normalized Distribution (Semi-Leptonic)
# ===================================================



plt.figure(figsize=(11, 12))
plt.subplots_adjust(left=0.15, right=0.95, bottom=0.12, top=0.95)

# ✅ Plot Signals
for signal_name, dsigma_data in signal_dsigma.items():
    pt_bins, dsigma = dsigma_data["pt_bins_lepton"]
    dsigma_norm = dsigma / (np.sum(dsigma) + 1e-12)
    plt.step(pt_bins, dsigma_norm, where="mid", alpha=0.7,
             label=f"Signal ($W^+W^-$) [{signal_name}]",
             color=signal_colors.get(signal_name, "black"), linewidth=3)
# ✅ Plot Backgrounds
background_styles = [
    (r"SM ($\gamma\gamma \to W^+W^-$)", background_dsigma_aa_ww, "blue", "-"),
    (r"$\gamma\gamma \to t\bar{t}$ $(\times 10^2)$", background_dsigma_aa_ttbar, "teal", "--"),
    (r"$\gamma\gamma \to \tau^+\tau^-$", background_dsigma_aa_tautau, "orange", "-"),
    (r"$\gamma\gamma \to \tau^+\tau^- (inel)$", background_dsigma_aa_tautau_inel, "brown", (0, (3, 1, 1, 1))),
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
    if np.sum(dsigma) > 0:
        dsigma_norm = dsigma / (np.sum(dsigma) + 1e-12)
        plt.step(pt_bins, dsigma_norm, where="mid", alpha=0.7,
                 label=label, color=color, linestyle=style, linewidth=3)
# ✅ Axis labels and title
plt.xlabel(r"$p_T^{\ell} \ \mathrm{[GeV]}$")
plt.ylabel(r"Normalized Distribution")
plt.title(r"Delphes simulation : $e^- p \to e^- W^+ W^- p \to e^- j j \ell \nu_{\ell} p$ : LHeC@1.2 TeV", fontsize=20)
#plt.yscale("log")
plt.ylim(1e-5, 0.4)
# ✅ Legend, Grid, and Save
plt.legend()
plt.grid(True, linestyle="--", alpha=0.6)
plt.tight_layout()
plt.savefig("/home/hamzeh-khanpour/Documents/GitHub/LHeC_Fast_Simulation/Pythia8_Delphes_semi_leptonic_allsignal_bkgs/normalized_distribution_lepton_pt_allFMsignal_allbkgs_all_decay_mode.pdf", dpi=600)

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
    (r"$\gamma\gamma \to \tau^+\tau^- (inel)$", background_dsigma_aa_tautau_inel, "brown", (0, (3, 1, 1, 1))),
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
plt.ylim(1e-5, 1e0)

# ✅ Legend and Layout
plt.legend()
plt.grid(True, linestyle="--", alpha=0.6)
plt.tight_layout()

# ✅ Save and Show
plt.savefig("/home/hamzeh-khanpour/Documents/GitHub/LHeC_Fast_Simulation/Pythia8_Delphes_semi_leptonic_allsignal_bkgs/differential_cross_section_jet_pt_allFMsignal_allbkgs_all_decay_mode.pdf", dpi=600)
plt.show()







# ===================================================
# ✅ Leading Jet \( p_T \) Normalized Distribution (Semi-Leptonic)
# ===================================================

plt.figure(figsize=(11, 12))
plt.subplots_adjust(left=0.15, right=0.95, bottom=0.12, top=0.95)

# ✅ Plot Signals
for signal_name, dsigma_data in signal_dsigma.items():
    pt_bins, dsigma = dsigma_data["pt_bins_jet"]
    dsigma_norm = dsigma / (np.sum(dsigma) + 1e-12)
    plt.step(pt_bins, dsigma_norm, where="mid", alpha=0.7,
             label=f"Signal ($W^+W^-$) [{signal_name}]",
             color=signal_colors.get(signal_name, "black"), linewidth=3)
# ✅ Plot Backgrounds
background_styles = [
    (r"SM ($\gamma\gamma \to W^+W^-$)", background_dsigma_aa_ww, "blue", "-"),
    (r"$\gamma\gamma \to t\bar{t}$ $(\times 10^2)$", background_dsigma_aa_ttbar, "teal", "--"),
    (r"$\gamma\gamma \to \tau^+\tau^-$", background_dsigma_aa_tautau, "orange", "-"),
    (r"$\gamma\gamma \to \tau^+\tau^- (inel)$", background_dsigma_aa_tautau_inel, "brown", (0, (3, 1, 1, 1))),
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
    if np.sum(dsigma) > 0:
        dsigma_norm = dsigma / (np.sum(dsigma) + 1e-12)
        plt.step(pt_bins, dsigma_norm, where="mid", alpha=0.7,
                 label=label, color=color, linestyle=style, linewidth=3)
# ✅ Axis Labels and Title
plt.xlabel(r"$p_T^{\mathrm{leading~jet}} \ \mathrm{[GeV]}$")
plt.ylabel(r"Normalized Distribution")
plt.title(r"Delphes simulation : $e^- p \to e^- W^+ W^- p \to e^- j j \ell \nu_{\ell} p$ : LHeC@1.2 \ TeV", fontsize=20)
#plt.yscale("log")
plt.ylim(1e-5, 0.20)
# ✅ Legend and Layout
plt.legend()
plt.grid(True, linestyle="--", alpha=0.6)
plt.tight_layout()
# ✅ Save and Show
plt.savefig("/home/hamzeh-khanpour/Documents/GitHub/LHeC_Fast_Simulation/Pythia8_Delphes_semi_leptonic_allsignal_bkgs/normalized_distribution_jet_pt_allFMsignal_allbkgs_all_decay_mode.pdf", dpi=600)
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
    (r"$\gamma\gamma \to \tau^+\tau^- (inel)$", background_dsigma_aa_tautau_inel, "brown", (0, (3, 1, 1, 1))),
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
# ✅ Lepton \( \eta \) Normalized Distribution (Semi-Leptonic)
# ===================================================

plt.figure(figsize=(11, 12))
plt.subplots_adjust(left=0.15, right=0.95, bottom=0.12, top=0.95)
# ✅ Plot Signals
for signal_name, dsigma_data in signal_dsigma.items():
    eta_bins, dsigma = dsigma_data["eta_bins_lepton"]
    dsigma_norm = dsigma / (np.sum(dsigma) + 1e-12)
    plt.step(eta_bins, dsigma_norm, where="mid", alpha=0.7,
             label=f"Signal ($W^+W^-$) [{signal_name}]",
             color=signal_colors.get(signal_name, "black"), linewidth=3)
# ✅ Plot Backgrounds
background_styles = [
    (r"SM ($\gamma\gamma \to W^+W^-$)", background_dsigma_aa_ww, "blue", "-"),
    (r"$\gamma\gamma \to t\bar{t}$ $(\times 10^2)$", background_dsigma_aa_ttbar, "teal", "--"),
    (r"$\gamma\gamma \to \tau^+\tau^-$", background_dsigma_aa_tautau, "orange", "-"),
    (r"$\gamma\gamma \to \tau^+\tau^- (inel)$", background_dsigma_aa_tautau_inel, "brown", (0, (3, 1, 1, 1))),
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
    if np.sum(dsigma) > 0:
        dsigma_norm = dsigma / (np.sum(dsigma) + 1e-12)
        plt.step(eta_bins, dsigma_norm, where="mid", alpha=0.7,
                 label=label, color=color, linestyle=style, linewidth=3)
# ✅ Axis Labels and Title
plt.xlabel(r"$\eta^{\ell}$")
plt.ylabel(r"Normalized Entries")
plt.title(r"Delphes simulation : $e^- p \to e^- W^+ W^- p \to e^- j j \ell \nu_{\ell} p$ : LHeC@1.2 TeV", fontsize=20)
#plt.yscale("log")
plt.ylim(1e-5, 0.2)
# ✅ Legend and Layout
plt.legend()
plt.grid(True, linestyle="--", alpha=0.6)
plt.tight_layout()
# ✅ Save and Show
plt.savefig("/home/hamzeh-khanpour/Documents/GitHub/LHeC_Fast_Simulation/Pythia8_Delphes_semi_leptonic_allsignal_bkgs/normalized_distribution_lepton_eta_allFMsignal_allbkgs_all_decay_mode.pdf", dpi=600)
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
    (r"$\gamma\gamma \to \tau^+\tau^- (inel)$", background_dsigma_aa_tautau_inel, "brown", (0, (3, 1, 1, 1))),
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
plt.ylim(1e-4, 100.0)

# ✅ Legend and Layout
plt.legend()
plt.grid(True, linestyle="--", alpha=0.6)
plt.tight_layout()

# ✅ Save and Show
plt.savefig("/home/hamzeh-khanpour/Documents/GitHub/LHeC_Fast_Simulation/Pythia8_Delphes_semi_leptonic_allsignal_bkgs/differential_cross_section_delta_r_allFMsignal_allbkgs_all_decay_mode.pdf", dpi=600)
plt.show()





# ===================================================
# ✅ Normalized ΔR(ℓ, jet) Distribution (Semi-Leptonic)
# ===================================================

plt.figure(figsize=(11, 12))
plt.subplots_adjust(left=0.15, right=0.95, bottom=0.12, top=0.95)

# ✅ Plot Normalized Signal Distributions
for signal_name, dsigma_data in signal_dsigma.items():
    delta_r_bins, dsigma = dsigma_data["delta_r_bins"]
    if np.sum(dsigma) > 0:
        normalized_dsigma = np.array(dsigma) / np.sum(dsigma)
        plt.step(delta_r_bins, normalized_dsigma, where="mid", alpha=0.7,
                 label=f"Signal ($W^+W^-$) [{signal_name}]",
                 color=signal_colors.get(signal_name, "black"), linewidth=3)

# ✅ Plot Normalized Backgrounds
background_styles = [
    (r"SM ($\gamma\gamma \to W^+W^-$)", background_dsigma_aa_ww, "blue", "-"),
    (r"$\gamma\gamma \to t\bar{t}$ $(\times 10^2)$", background_dsigma_aa_ttbar, "teal", "--"),
    (r"$\gamma\gamma \to \tau^+\tau^-$", background_dsigma_aa_tautau, "orange", "-"),
    (r"$\gamma\gamma \to \tau^+\tau^- (inel)$", background_dsigma_aa_tautau_inel, "brown", (0, (3, 1, 1, 1))),
    (r"Inclusive $t\bar{t}$", background_dsigma_inclusive_ttbar, "orchid", (0, (3, 1, 1, 1))),
    (r"Single Top", background_dsigma_single_top, "teal", "--"),
    (r"W Production", background_dsigma_w_production, "darkgreen", "-."),
    (r"Z Production", background_dsigma_z_production, "darkred", ":"),
    (r"$WWj$", background_dsigma_wwj_production, "royalblue", "--"),
    (r"$ZZj$ $(\times 10^2)$", background_dsigma_zzj_production, "mediumvioletred", "-."),
    (r"$WZj$", background_dsigma_wzj_production, "slategray", ":")
]

for label, dsigma_dict, color, style in background_styles:
    delta_r_bins, dsigma = dsigma_dict["delta_r_bins"]
    if np.sum(dsigma) > 0:
        normalized_dsigma = np.array(dsigma) / np.sum(dsigma)
        plt.step(delta_r_bins, normalized_dsigma, where="mid", alpha=0.7,
                 label=label, color=color, linestyle=style, linewidth=3)

# ✅ Axis Labels and Title
plt.xlabel(r"$\Delta R(\ell, \mathrm{leading~jet})$")
plt.ylabel(r"Normalized Entries")
plt.title(r"Delphes simulation : $e^- p \to e^- W^+ W^- p \to e^- j j \ell \nu_{\ell} p$ : LHeC@1.2 TeV", fontsize=20)
plt.ylim(1e-5, 0.3)
#plt.yscale("log")
plt.legend()
plt.grid(True, linestyle="--", alpha=0.6)
plt.tight_layout()
plt.savefig("/home/hamzeh-khanpour/Documents/GitHub/LHeC_Fast_Simulation/Pythia8_Delphes_semi_leptonic_allsignal_bkgs/normalized_delta_r_allFMsignal_allbkgs_all_decay_mode.pdf", dpi=600)
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
    (r"$\gamma\gamma \to \tau^+\tau^- (inel)$", background_dsigma_aa_tautau_inel, "brown", (0, (3, 1, 1, 1))),
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
plt.ylim(1e-5, 10.0)

# ✅ Legend and Layout
plt.legend()
plt.grid(True, linestyle="--", alpha=0.6)
plt.tight_layout()

# ✅ Save and Show
plt.savefig("/home/hamzeh-khanpour/Documents/GitHub/LHeC_Fast_Simulation/Pythia8_Delphes_semi_leptonic_allsignal_bkgs/differential_cross_section_met_allFMsignal_allbkgs_all_decay_mode.pdf", dpi=600)
plt.show()






# ===================================================
# ✅ Normalized MET Distribution (Semi-Leptonic)
# ===================================================

plt.figure(figsize=(11, 12))
plt.subplots_adjust(left=0.15, right=0.95, bottom=0.12, top=0.95)

# ✅ Normalized Signal Distributions
for signal_name, dsigma_data in signal_dsigma.items():
    met_bins, dsigma = dsigma_data["met_bins"]
    if np.sum(dsigma) > 0:
        normalized_dsigma = np.array(dsigma) / np.sum(dsigma)
        plt.step(met_bins, normalized_dsigma, where="mid", alpha=0.7,
                 label=f"Signal ($W^+W^-$) [{signal_name}]",
                 color=signal_colors.get(signal_name, "black"), linewidth=3)

# ✅ Normalized Backgrounds
background_styles = [
    (r"SM ($\gamma\gamma \to W^+W^-$)", background_dsigma_aa_ww, "blue", "-"),
    (r"$\gamma\gamma \to t\bar{t}$ $(\times 10^2)$", background_dsigma_aa_ttbar, "teal", "--"),
    (r"$\gamma\gamma \to \tau^+\tau^-$", background_dsigma_aa_tautau, "orange", "-"),
    (r"$\gamma\gamma \to \tau^+\tau^- (inel)$", background_dsigma_aa_tautau_inel, "brown", (0, (3, 1, 1, 1))),
    (r"Inclusive $t\bar{t}$", background_dsigma_inclusive_ttbar, "orchid", (0, (3, 1, 1, 1))),
    (r"Single Top", background_dsigma_single_top, "teal", "--"),
    (r"W Production", background_dsigma_w_production, "darkgreen", "-."),
    (r"Z Production", background_dsigma_z_production, "darkred", ":"),
    (r"$WWj$", background_dsigma_wwj_production, "royalblue", "--"),
    (r"$ZZj$ $(\times 10^2)$", background_dsigma_zzj_production, "mediumvioletred", "-."),
    (r"$WZj$", background_dsigma_wzj_production, "slategray", ":")
]

for label, dsigma_dict, color, style in background_styles:
    met_bins, dsigma = dsigma_dict["met_bins"]
    if np.sum(dsigma) > 0:
        normalized_dsigma = np.array(dsigma) / np.sum(dsigma)
        plt.step(met_bins, normalized_dsigma, where="mid", alpha=0.7,
                 label=label, color=color, linestyle=style, linewidth=3)

# ✅ Axis Labels and Title
plt.xlabel(r"$\mathrm{MET} \ \mathrm{[GeV]}$")
plt.ylabel(r"Normalized Entries")
plt.title(r"Delphes simulation : $e^- p \to e^- W^+ W^- p \to e^- j j \ell \nu_{\ell} p$ : LHeC@1.2 TeV", fontsize=20)
#plt.yscale("log")
plt.ylim(1e-5, 0.3)

# ✅ Legend and Save
plt.legend()
plt.grid(True, linestyle="--", alpha=0.6)
plt.tight_layout()
plt.savefig("/home/hamzeh-khanpour/Documents/GitHub/LHeC_Fast_Simulation/Pythia8_Delphes_semi_leptonic_allsignal_bkgs/normalized_met_allFMsignal_allbkgs_all_decay_mode.pdf", dpi=600)
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
    (r"$\gamma\gamma \to \tau^+\tau^- (inel)$", background_dsigma_aa_tautau_inel, "brown", (0, (3, 1, 1, 1))),
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
    (r"$\gamma\gamma \to \tau^+\tau^- (inel)$", background_dsigma_aa_tautau_inel, "brown", (0, (3, 1, 1, 1))),
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
plt.ylim(0.0, 0.3)

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
    (r"$\gamma\gamma \to \tau^+\tau^- (inel)$", background_dsigma_aa_tautau_inel, "brown", (0, (3, 1, 1, 1))),
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
    (r"$\gamma\gamma \to \tau^+\tau^- (inel)$", background_dsigma_aa_tautau_inel, "brown", (0, (3, 1, 1, 1))),
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
    (r"$\gamma\gamma \to \tau^+\tau^- (inel)$", background_dsigma_aa_tautau_inel, "brown", (0, (3, 1, 1, 1))),
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








# ===================================================
# ✅ Normalized Δη_{jj} Distribution (Semi-Leptonic)
# ===================================================

plt.figure(figsize=(11, 12))
plt.subplots_adjust(left=0.15, right=0.95, bottom=0.12, top=0.95)

# ✅ Plot Normalized Signal Distributions
for signal_name, dsigma_data in signal_dsigma.items():
    delta_eta_jj_bins, dsigma = dsigma_data["delta_eta_jj_bins"]
    if np.sum(dsigma) > 0:
        normalized_dsigma = dsigma / np.sum(dsigma)
        plt.step(delta_eta_jj_bins, normalized_dsigma, where="mid", alpha=0.7,
                 label=f"Signal ($W^+W^-$) [{signal_name}]",
                 color=signal_colors.get(signal_name, "black"), linewidth=3)

# ✅ Plot Normalized Backgrounds
background_styles = [
    (r"SM ($\gamma\gamma \to W^+W^-$)", background_dsigma_aa_ww, "blue", "-"),
    (r"$\gamma\gamma \to t\bar{t}$ $(\times 10^2)$", background_dsigma_aa_ttbar, "teal", "--"),
    (r"$\gamma\gamma \to \tau^+\tau^-$", background_dsigma_aa_tautau, "orange", "-"),
    (r"$\gamma\gamma \to \tau^+\tau^- (inel)$", background_dsigma_aa_tautau_inel, "brown", (0, (3, 1, 1, 1))),
    (r"Inclusive $t\bar{t}$", background_dsigma_inclusive_ttbar, "orchid", (0, (3, 1, 1, 1))),
    (r"Single Top", background_dsigma_single_top, "teal", "--"),
    (r"W Production", background_dsigma_w_production, "darkgreen", "-."),
    (r"Z Production", background_dsigma_z_production, "darkred", ":"),
    (r"$WWj$", background_dsigma_wwj_production, "royalblue", (0, (5, 2))),
    (r"$ZZj$ $(\times 10^2)$", background_dsigma_zzj_production, "crimson", (0, (1, 1))),
    (r"$WZj$", background_dsigma_wzj_production, "slateblue", (0, (3, 2, 1, 2)))
]

for label, dsigma_dict, color, style in background_styles:
    delta_eta_jj_bins, dsigma = dsigma_dict["delta_eta_jj_bins"]
    if np.sum(dsigma) > 0:
        normalized_dsigma = dsigma / np.sum(dsigma)
        plt.step(delta_eta_jj_bins, normalized_dsigma, where="mid", alpha=0.7,
                 label=label, color=color, linestyle=style, linewidth=3)

# ✅ Axis Labels and Title
plt.xlabel(r"$\Delta\eta(j_1, j_2)$")
plt.ylabel("Normalized Entries")
plt.title(r"Delphes simulation : $e^- p \to e^- W^+ W^- p \to e^- j j \ell \nu_{\ell} p$ : LHeC@1.2 TeV", fontsize=20)
plt.ylim(0.0, 0.20)

# ✅ Legend and Layout
plt.legend()
plt.grid(True, linestyle="--", alpha=0.6)
plt.tight_layout()

# ✅ Save and Show
plt.savefig("/home/hamzeh-khanpour/Documents/GitHub/LHeC_Fast_Simulation/Pythia8_Delphes_semi_leptonic_allsignal_bkgs/normalized_delta_eta_jj_allFMsignal_allbkgs_all_decay_mode.pdf", dpi=600)
plt.show()







# ===================================================
# ✅ Leptonic W \( p_T \) Differential Cross-Section (Semi-Leptonic)
# ===================================================

plt.figure(figsize=(11, 12))  # Create a new figure for W_lep pT
plt.subplots_adjust(left=0.15, right=0.95, bottom=0.12, top=0.95)

# ✅ Plot Signals
for signal_name, dsigma_data in signal_dsigma.items():
    pt_bins, dsigma = dsigma_data["pt_w_leptonic_bins"]
    plt.step(pt_bins, dsigma, where="mid", alpha=0.7,
             label=f"Signal ($W^+W^-$) [{signal_name}]",
             color=signal_colors.get(signal_name, "black"), linewidth=3)

# ✅ Plot Backgrounds
background_styles = [
    (r"SM ($\gamma\gamma \to W^+W^-$)", background_dsigma_aa_ww, "blue", "-"),
    (r"$\gamma\gamma \to t\bar{t}$ $(\times 10^2)$", background_dsigma_aa_ttbar, "teal", "--"),
    (r"$\gamma\gamma \to \tau^+\tau^-$", background_dsigma_aa_tautau, "orange", "-"),
    (r"$\gamma\gamma \to \tau^+\tau^- (inel)$", background_dsigma_aa_tautau_inel, "brown", (0, (3, 1, 1, 1))),
    (r"Inclusive $t\bar{t}$", background_dsigma_inclusive_ttbar, "orchid", (0, (3, 1, 1, 1))),
    (r"Single Top", background_dsigma_single_top, "teal", "--"),
    (r"W Production", background_dsigma_w_production, "darkgreen", "-."),
    (r"Z Production", background_dsigma_z_production, "darkred", ":"),
    (r"$WWj$", background_dsigma_wwj_production, "royalblue", "--"),
    (r"$ZZj$ $(\times 10^2)$", background_dsigma_zzj_production, "mediumvioletred", "-."),
    (r"$WZj$", background_dsigma_wzj_production, "slategray", ":")
]

for label, dsigma_dict, color, style in background_styles:
    pt_bins, dsigma = dsigma_dict["pt_w_leptonic_bins"]
    if sum(dsigma) > 0:
        plt.step(pt_bins, dsigma, where="mid", alpha=0.7,
                 label=label, color=color, linestyle=style, linewidth=3)

# ✅ Axis Labels and Title
plt.xlabel(r"$p_T^{W^{\mathrm{lep}}} \ \mathrm{[GeV]}$")
plt.ylabel(r"$\frac{d\sigma}{dp_T^{W^{\mathrm{lep}}}} \ \mathrm{[pb/GeV]}$")
plt.title(r"Delphes simulation : $e^- p \to e^- W^+ W^- p \to e^- j j \ell \nu_{\ell} p$ : LHeC@1.2 TeV", fontsize=20)
plt.yscale("log")
plt.ylim(1e-5, 1e0)

# ✅ Legend, Grid, Save
plt.legend()
plt.grid(True, linestyle="--", alpha=0.6)
plt.tight_layout()
plt.savefig("/home/hamzeh-khanpour/Documents/GitHub/LHeC_Fast_Simulation/Pythia8_Delphes_semi_leptonic_allsignal_bkgs/differential_cross_section_pt_w_leptonic_allFMsignal_allbkgs.pdf", dpi=600)

plt.show()







# ===================================================
# ✅ Leptonic W \( p_T \) Differential Cross-Section (Semi-Leptonic)
# ===================================================

plt.figure(figsize=(11, 12))  # Create a new figure for W_lep pT
plt.subplots_adjust(left=0.15, right=0.95, bottom=0.12, top=0.95)

# ✅ Plot Signals
for signal_name, dsigma_data in signal_dsigma.items():
    pt_bins, dsigma = dsigma_data["pt_w_leptonic_bins"]
    total = sum(dsigma)
    if total > 0:
        dsigma_norm = [x / total for x in dsigma]
        plt.step(pt_bins, dsigma_norm, where="mid", alpha=0.7,
                 label=f"Signal ($W^+W^-$) [{signal_name}]",
                 color=signal_colors.get(signal_name, "black"), linewidth=3)

# ✅ Plot Backgrounds
background_styles = [
    (r"SM ($\gamma\gamma \to W^+W^-$)", background_dsigma_aa_ww, "blue", "-"),
    (r"$\gamma\gamma \to t\bar{t}$ $(\times 10^2)$", background_dsigma_aa_ttbar, "teal", "--"),
    (r"$\gamma\gamma \to \tau^+\tau^-$", background_dsigma_aa_tautau, "orange", "-"),
    (r"$\gamma\gamma \to \tau^+\tau^- (inel)$", background_dsigma_aa_tautau_inel, "brown", (0, (3, 1, 1, 1))),
    (r"Inclusive $t\bar{t}$", background_dsigma_inclusive_ttbar, "orchid", (0, (3, 1, 1, 1))),
    (r"Single Top", background_dsigma_single_top, "teal", "--"),
    (r"W Production", background_dsigma_w_production, "darkgreen", "-."),
    (r"Z Production", background_dsigma_z_production, "darkred", ":"),
    (r"$WWj$", background_dsigma_wwj_production, "royalblue", "--"),
    (r"$ZZj$ $(\times 10^2)$", background_dsigma_zzj_production, "mediumvioletred", "-."),
    (r"$WZj$", background_dsigma_wzj_production, "slategray", ":")
]

for label, dsigma_dict, color, style in background_styles:
    pt_bins, dsigma = dsigma_dict["pt_w_leptonic_bins"]
    total = sum(dsigma)
    if total > 0:
        dsigma_norm = [x / total for x in dsigma]
        plt.step(pt_bins, dsigma_norm, where="mid", alpha=0.7,
                 label=label, color=color, linestyle=style, linewidth=3)

# ✅ Axis Labels and Title
plt.xlabel(r"$p_T^{W^{\mathrm{lep}}} \ \mathrm{[GeV]}$")
plt.ylabel(r"Normalized Entries")
plt.title(r"Delphes simulation : $e^- p \to e^- W^+ W^- p \to e^- j j \ell \nu_{\ell} p$ : LHeC@1.2 TeV", fontsize=20)
#plt.yscale("log")
plt.ylim(1e-5, 0.3)

# ✅ Legend, Grid, Save
plt.legend()
plt.grid(True, linestyle="--", alpha=0.6)
plt.tight_layout()
plt.savefig("/home/hamzeh-khanpour/Documents/GitHub/LHeC_Fast_Simulation/Pythia8_Delphes_semi_leptonic_allsignal_bkgs/normalized_differential_cross_section_pt_w_leptonic_allFMsignal_allbkgs.pdf", dpi=600)

plt.show()










# ===================================================
# ✅ Hadronic W \( p_T \) Differential Cross-Section (Semi-Leptonic)
# ===================================================

plt.figure(figsize=(11, 12))  # New figure for hadronic W pT
plt.subplots_adjust(left=0.15, right=0.95, bottom=0.12, top=0.95)

# ✅ Plot Signals
for signal_name, dsigma_data in signal_dsigma.items():
    pt_bins, dsigma = dsigma_data["pt_w_hadronic_bins"]
    plt.step(pt_bins, dsigma, where="mid", alpha=0.7,
             label=f"Signal ($W^+W^-$) [{signal_name}]",
             color=signal_colors.get(signal_name, "black"), linewidth=3)

# ✅ Plot Backgrounds
background_styles = [
    (r"SM ($\gamma\gamma \to W^+W^-$)", background_dsigma_aa_ww, "blue", "-"),
    (r"$\gamma\gamma \to t\bar{t}$ $(\times 10^2)$", background_dsigma_aa_ttbar, "teal", "--"),
    (r"$\gamma\gamma \to \tau^+\tau^-$", background_dsigma_aa_tautau, "orange", "-"),
    (r"$\gamma\gamma \to \tau^+\tau^- (inel)$", background_dsigma_aa_tautau_inel, "brown", (0, (3, 1, 1, 1))),
    (r"Inclusive $t\bar{t}$", background_dsigma_inclusive_ttbar, "orchid", (0, (3, 1, 1, 1))),
    (r"Single Top", background_dsigma_single_top, "teal", "--"),
    (r"W Production", background_dsigma_w_production, "darkgreen", "-."),
    (r"Z Production", background_dsigma_z_production, "darkred", ":"),
    (r"$WWj$", background_dsigma_wwj_production, "royalblue", "--"),
    (r"$ZZj$ $(\times 10^2)$", background_dsigma_zzj_production, "mediumvioletred", "-."),
    (r"$WZj$", background_dsigma_wzj_production, "slategray", ":")
]

for label, dsigma_dict, color, style in background_styles:
    pt_bins, dsigma = dsigma_dict["pt_w_hadronic_bins"]
    if sum(dsigma) > 0:
        plt.step(pt_bins, dsigma, where="mid", alpha=0.7,
                 label=label, color=color, linestyle=style, linewidth=3)

# ✅ Axis Labels and Title
plt.xlabel(r"$p_T^{W^{\mathrm{had}}} \ \mathrm{[GeV]}$")
plt.ylabel(r"$\frac{d\sigma}{dp_T^{W^{\mathrm{had}}}} \ \mathrm{[pb/GeV]}$")
plt.title(r"Delphes simulation : $e^- p \to e^- W^+ W^- p \to e^- j j \ell \nu_{\ell} p$ : LHeC@1.2 TeV", fontsize=20)
plt.yscale("log")
plt.ylim(1e-5, 1e0)

# ✅ Legend, Grid, Save
plt.legend()
plt.grid(True, linestyle="--", alpha=0.6)
plt.tight_layout()
plt.savefig("/home/hamzeh-khanpour/Documents/GitHub/LHeC_Fast_Simulation/Pythia8_Delphes_semi_leptonic_allsignal_bkgs/differential_cross_section_pt_w_hadronic_allFMsignal_allbkgs.pdf", dpi=600)

plt.show()







# ===================================================
# ✅ Hadronic W \( p_T \) Differential Cross-Section (Semi-Leptonic)
# ===================================================

plt.figure(figsize=(11, 12))  # New figure for hadronic W pT
plt.subplots_adjust(left=0.15, right=0.95, bottom=0.12, top=0.95)

# ✅ Plot Signals
for signal_name, dsigma_data in signal_dsigma.items():
    pt_bins, dsigma = dsigma_data["pt_w_hadronic_bins"]
    total = sum(dsigma)
    if total > 0:
        dsigma_norm = [x / total for x in dsigma]
        plt.step(pt_bins, dsigma_norm, where="mid", alpha=0.7,
                 label=f"Signal ($W^+W^-$) [{signal_name}]",
                 color=signal_colors.get(signal_name, "black"), linewidth=3)

# ✅ Plot Backgrounds
background_styles = [
    (r"SM ($\gamma\gamma \to W^+W^-$)", background_dsigma_aa_ww, "blue", "-"),
    (r"$\gamma\gamma \to t\bar{t}$ $(\times 10^2)$", background_dsigma_aa_ttbar, "teal", "--"),
    (r"$\gamma\gamma \to \tau^+\tau^-$", background_dsigma_aa_tautau, "orange", "-"),
    (r"$\gamma\gamma \to \tau^+\tau^- (inel)$", background_dsigma_aa_tautau_inel, "brown", (0, (3, 1, 1, 1))),
    (r"Inclusive $t\bar{t}$", background_dsigma_inclusive_ttbar, "orchid", (0, (3, 1, 1, 1))),
    (r"Single Top", background_dsigma_single_top, "teal", "--"),
    (r"W Production", background_dsigma_w_production, "darkgreen", "-."),
    (r"Z Production", background_dsigma_z_production, "darkred", ":"),
    (r"$WWj$", background_dsigma_wwj_production, "royalblue", "--"),
    (r"$ZZj$ $(\times 10^2)$", background_dsigma_zzj_production, "mediumvioletred", "-."),
    (r"$WZj$", background_dsigma_wzj_production, "slategray", ":")
]

for label, dsigma_dict, color, style in background_styles:
    pt_bins, dsigma = dsigma_dict["pt_w_hadronic_bins"]
    total = sum(dsigma)
    if total > 0:
        dsigma_norm = [x / total for x in dsigma]
        plt.step(pt_bins, dsigma_norm, where="mid", alpha=0.7,
                 label=label, color=color, linestyle=style, linewidth=3)

# ✅ Axis Labels and Title
plt.xlabel(r"$p_T^{W^{\mathrm{had}}} \ \mathrm{[GeV]}$")
plt.ylabel(r"Normalized Entries")
plt.title(r"Delphes simulation : $e^- p \to e^- W^+ W^- p \to e^- j j \ell \nu_{\ell} p$ : LHeC@1.2 TeV", fontsize=20)
#plt.yscale("log")
plt.ylim(1e-5, 0.3)

# ✅ Legend, Grid, Save
plt.legend()
plt.grid(True, linestyle="--", alpha=0.6)
plt.tight_layout()
plt.savefig("/home/hamzeh-khanpour/Documents/GitHub/LHeC_Fast_Simulation/Pythia8_Delphes_semi_leptonic_allsignal_bkgs/normalized_differential_cross_section_pt_w_hadronic_allFMsignal_allbkgs.pdf", dpi=600)

plt.show()









# ===================================================
# ✅ Δϕ(ℓ, MET) Differential Cross-Section (Semi-Leptonic)
# ===================================================

plt.figure(figsize=(11, 12))  # Create a new figure for Δϕ(l, MET)
plt.subplots_adjust(left=0.15, right=0.95, bottom=0.12, top=0.95)

# ✅ Plot Signals
for signal_name, dsigma_data in signal_dsigma.items():
    delta_phi_bins, dsigma = dsigma_data["delta_phi_lep_met_bins"]
    plt.step(delta_phi_bins, dsigma, where="mid", alpha=0.7,
             label=f"Signal ($W^+W^-$) [{signal_name}]",
             color=signal_colors.get(signal_name, "black"), linewidth=3)

# ✅ Plot Backgrounds
background_styles = [
    (r"SM ($\gamma\gamma \to W^+W^-$)", background_dsigma_aa_ww, "blue", "-"),
    (r"$\gamma\gamma \to t\bar{t}$ $(\times 10^2)$", background_dsigma_aa_ttbar, "teal", "--"),
    (r"$\gamma\gamma \to \tau^+\tau^-$", background_dsigma_aa_tautau, "orange", "-"),
    (r"$\gamma\gamma \to \tau^+\tau^- (inel)$", background_dsigma_aa_tautau_inel, "brown", (0, (3, 1, 1, 1))),
    (r"Inclusive $t\bar{t}$", background_dsigma_inclusive_ttbar, "orchid", (0, (3, 1, 1, 1))),
    (r"Single Top", background_dsigma_single_top, "teal", "--"),
    (r"W Production", background_dsigma_w_production, "darkgreen", "-."),
    (r"Z Production", background_dsigma_z_production, "darkred", ":"),
    (r"$WWj$", background_dsigma_wwj_production, "royalblue", "--"),
    (r"$ZZj$ $(\times 10^2)$", background_dsigma_zzj_production, "mediumvioletred", "-."),
    (r"$WZj$", background_dsigma_wzj_production, "slategray", ":")
]

for label, dsigma_dict, color, style in background_styles:
    delta_phi_bins, dsigma = dsigma_dict["delta_phi_lep_met_bins"]
    if sum(dsigma) > 0:
        plt.step(delta_phi_bins, dsigma, where="mid", alpha=0.7,
                 label=label, color=color, linestyle=style, linewidth=3)

# ✅ Labels and Title
plt.xlabel(r"$\Delta\phi(\ell, \mathrm{MET})$")
plt.ylabel(r"$\frac{d\sigma}{d\Delta\phi(\ell, \mathrm{MET})} \ \mathrm{[pb]}$")
plt.title(r"Delphes simulation : $e^- p \to e^- W^+ W^- p \to e^- j j \ell \nu_{\ell} p$ : LHeC@1.2 TeV", fontsize=20)
plt.yscale("log")
plt.ylim(1e-5, 1e1)

# ✅ Legend, Grid, Save
plt.legend()
plt.grid(True, linestyle="--", alpha=0.6)
plt.tight_layout()
plt.savefig("/home/hamzeh-khanpour/Documents/GitHub/LHeC_Fast_Simulation/Pythia8_Delphes_semi_leptonic_allsignal_bkgs/dsigma_delta_phi_lep_met_allFMsignal_allbkgs.pdf", dpi=600)

plt.show()






# ===================================================
# ✅ Δϕ(ℓ, MET) Differential Cross-Section (Semi-Leptonic)
# ===================================================

plt.figure(figsize=(11, 12))  # Create a new figure for Δϕ(l, MET)
plt.subplots_adjust(left=0.15, right=0.95, bottom=0.12, top=0.95)

# ✅ Plot Signals
for signal_name, dsigma_data in signal_dsigma.items():
    delta_phi_bins, dsigma = dsigma_data["delta_phi_lep_met_bins"]
    total = sum(dsigma)
    if total > 0:
        dsigma_norm = [x / total for x in dsigma]
        plt.step(delta_phi_bins, dsigma_norm, where="mid", alpha=0.7,
                 label=f"Signal ($W^+W^-$) [{signal_name}]",
                 color=signal_colors.get(signal_name, "black"), linewidth=3)

# ✅ Plot Backgrounds
background_styles = [
    (r"SM ($\gamma\gamma \to W^+W^-$)", background_dsigma_aa_ww, "blue", "-"),
    (r"$\gamma\gamma \to t\bar{t}$ $(\times 10^2)$", background_dsigma_aa_ttbar, "teal", "--"),
    (r"$\gamma\gamma \to \tau^+\tau^-$", background_dsigma_aa_tautau, "orange", "-"),
    (r"$\gamma\gamma \to \tau^+\tau^- (inel)$", background_dsigma_aa_tautau_inel, "brown", (0, (3, 1, 1, 1))),
    (r"Inclusive $t\bar{t}$", background_dsigma_inclusive_ttbar, "orchid", (0, (3, 1, 1, 1))),
    (r"Single Top", background_dsigma_single_top, "teal", "--"),
    (r"W Production", background_dsigma_w_production, "darkgreen", "-."),
    (r"Z Production", background_dsigma_z_production, "darkred", ":"),
    (r"$WWj$", background_dsigma_wwj_production, "royalblue", "--"),
    (r"$ZZj$ $(\times 10^2)$", background_dsigma_zzj_production, "mediumvioletred", "-."),
    (r"$WZj$", background_dsigma_wzj_production, "slategray", ":")
]

for label, dsigma_dict, color, style in background_styles:
    delta_phi_bins, dsigma = dsigma_dict["delta_phi_lep_met_bins"]
    total = sum(dsigma)
    if total > 0:
        dsigma_norm = [x / total for x in dsigma]
        plt.step(delta_phi_bins, dsigma_norm, where="mid", alpha=0.7,
                 label=label, color=color, linestyle=style, linewidth=3)

# ✅ Labels and Title
plt.xlabel(r"$\Delta\phi(\ell, \mathrm{MET})$")
plt.ylabel(r"Normalized Entries")
plt.title(r"Delphes simulation : $e^- p \to e^- W^+ W^- p \to e^- j j \ell \nu_{\ell} p$ : LHeC@1.2 TeV", fontsize=20)
#plt.yscale("log")
plt.ylim(1e-5, 0.2)

# ✅ Legend, Grid, Save
plt.legend()
plt.grid(True, linestyle="--", alpha=0.6)
plt.tight_layout()
plt.savefig("/home/hamzeh-khanpour/Documents/GitHub/LHeC_Fast_Simulation/Pythia8_Delphes_semi_leptonic_allsignal_bkgs/normalized_dsigma_delta_phi_lep_met_allFMsignal_allbkgs.pdf", dpi=600)

plt.show()










# ===================================================
# ✅ Transverse Mass of Leptonic W (M_T) (Semi-Leptonic)
# ===================================================

plt.figure(figsize=(11, 12))  # New figure for MT(W_lep)
plt.subplots_adjust(left=0.15, right=0.95, bottom=0.12, top=0.95)

# ✅ Plot Signals
for signal_name, dsigma_data in signal_dsigma.items():
    mt_bins, dsigma = dsigma_data["mt_w_leptonic_bins"]
    plt.step(mt_bins, dsigma, where="mid", alpha=0.7,
             label=f"Signal ($W^+W^-$) [{signal_name}]",
             color=signal_colors.get(signal_name, "black"), linewidth=3)

# ✅ Plot Backgrounds
background_styles = [
    (r"SM ($\gamma\gamma \to W^+W^-$)", background_dsigma_aa_ww, "blue", "-"),
    (r"$\gamma\gamma \to t\bar{t}$ $(\times 10^2)$", background_dsigma_aa_ttbar, "teal", "--"),
    (r"$\gamma\gamma \to \tau^+\tau^-$", background_dsigma_aa_tautau, "orange", "-"),
    (r"$\gamma\gamma \to \tau^+\tau^- (inel)$", background_dsigma_aa_tautau_inel, "brown", (0, (3, 1, 1, 1))),
    (r"Inclusive $t\bar{t}$", background_dsigma_inclusive_ttbar, "orchid", (0, (3, 1, 1, 1))),
    (r"Single Top", background_dsigma_single_top, "teal", "--"),
    (r"W Production", background_dsigma_w_production, "darkgreen", "-."),
    (r"Z Production", background_dsigma_z_production, "darkred", ":"),
    (r"$WWj$", background_dsigma_wwj_production, "royalblue", "--"),
    (r"$ZZj$ $(\times 10^2)$", background_dsigma_zzj_production, "mediumvioletred", "-."),
    (r"$WZj$", background_dsigma_wzj_production, "slategray", ":")
]

for label, dsigma_dict, color, style in background_styles:
    mt_bins, dsigma = dsigma_dict["mt_w_leptonic_bins"]
    if sum(dsigma) > 0:
        plt.step(mt_bins, dsigma, where="mid", alpha=0.7,
                 label=label, color=color, linestyle=style, linewidth=3)

# ✅ Axis labels and layout
plt.xlabel(r"$M_T(W^\mathrm{lep}) \ \mathrm{[GeV]}$")
plt.ylabel(r"$\frac{d\sigma}{dM_T} \ \mathrm{[pb/GeV]}$")
plt.title(r"Delphes simulation : $e^- p \to e^- W^+ W^- p \to e^- j j \ell \nu_{\ell} p$ : LHeC@1.2 TeV", fontsize=20)
plt.yscale("log")
plt.ylim(1e-5, 1e1)

# ✅ Legend and Save
plt.legend()
plt.grid(True, linestyle="--", alpha=0.6)
plt.tight_layout()
plt.savefig("/home/hamzeh-khanpour/Documents/GitHub/LHeC_Fast_Simulation/Pythia8_Delphes_semi_leptonic_allsignal_bkgs/dsigma_mt_w_leptonic_allFMsignal_allbkgs.pdf", dpi=600)

plt.show()






# ===================================================
# ✅ Transverse Mass of Leptonic W (M_T) (Semi-Leptonic)
# ===================================================

plt.figure(figsize=(11, 12))  # New figure for MT(W_lep)
plt.subplots_adjust(left=0.15, right=0.95, bottom=0.12, top=0.95)

# ✅ Plot Signals
for signal_name, dsigma_data in signal_dsigma.items():
    mt_bins, dsigma = dsigma_data["mt_w_leptonic_bins"]
    total = sum(dsigma)
    if total > 0:
        dsigma_norm = [x / total for x in dsigma]
        plt.step(mt_bins, dsigma_norm, where="mid", alpha=0.7,
                 label=f"Signal ($W^+W^-$) [{signal_name}]",
                 color=signal_colors.get(signal_name, "black"), linewidth=3)

# ✅ Plot Backgrounds
background_styles = [
    (r"SM ($\gamma\gamma \to W^+W^-$)", background_dsigma_aa_ww, "blue", "-"),
    (r"$\gamma\gamma \to t\bar{t}$ $(\times 10^2)$", background_dsigma_aa_ttbar, "teal", "--"),
    (r"$\gamma\gamma \to \tau^+\tau^-$", background_dsigma_aa_tautau, "orange", "-"),
    (r"$\gamma\gamma \to \tau^+\tau^- (inel)$", background_dsigma_aa_tautau_inel, "brown", (0, (3, 1, 1, 1))),
    (r"Inclusive $t\bar{t}$", background_dsigma_inclusive_ttbar, "orchid", (0, (3, 1, 1, 1))),
    (r"Single Top", background_dsigma_single_top, "teal", "--"),
    (r"W Production", background_dsigma_w_production, "darkgreen", "-."),
    (r"Z Production", background_dsigma_z_production, "darkred", ":"),
    (r"$WWj$", background_dsigma_wwj_production, "royalblue", "--"),
    (r"$ZZj$ $(\times 10^2)$", background_dsigma_zzj_production, "mediumvioletred", "-."),
    (r"$WZj$", background_dsigma_wzj_production, "slategray", ":")
]

for label, dsigma_dict, color, style in background_styles:
    mt_bins, dsigma = dsigma_dict["mt_w_leptonic_bins"]
    total = sum(dsigma)
    if total > 0:
        dsigma_norm = [x / total for x in dsigma]
        plt.step(mt_bins, dsigma_norm, where="mid", alpha=0.7,
                 label=label, color=color, linestyle=style, linewidth=3)

# ✅ Axis labels and layout
plt.xlabel(r"$M_T(W^\mathrm{lep}) \ \mathrm{[GeV]}$")
plt.ylabel(r"Normalized $\frac{1}{\sigma} \frac{d\sigma}{dM_T}$")
plt.title(r"Delphes simulation : $e^- p \to e^- W^+ W^- p \to e^- j j \ell \nu_{\ell} p$ : LHeC@1.2 TeV", fontsize=20)
#plt.yscale("log")
plt.ylim(1e-5, 0.3)

# ✅ Legend and Save
plt.legend()
plt.grid(True, linestyle="--", alpha=0.6)
plt.tight_layout()
plt.savefig("/home/hamzeh-khanpour/Documents/GitHub/LHeC_Fast_Simulation/Pythia8_Delphes_semi_leptonic_allsignal_bkgs/normalized_dsigma_mt_w_leptonic_allFMsignal_allbkgs.pdf", dpi=600)

plt.show()










# ===================================================
# ✅ Total Scalar Transverse Energy \( H_T \) (Semi-Leptonic)
# ===================================================

plt.figure(figsize=(11, 12))  # New figure for H_T
plt.subplots_adjust(left=0.15, right=0.95, bottom=0.12, top=0.95)

# ✅ Plot Signals
for signal_name, dsigma_data in signal_dsigma.items():
    ht_bins, dsigma = dsigma_data["ht_total_bins"]
    plt.step(ht_bins, dsigma, where="mid", alpha=0.7,
             label=f"Signal ($W^+W^-$) [{signal_name}]",
             color=signal_colors.get(signal_name, "black"), linewidth=3)

# ✅ Plot Backgrounds
background_styles = [
    (r"SM ($\gamma\gamma \to W^+W^-$)", background_dsigma_aa_ww, "blue", "-"),
    (r"$\gamma\gamma \to t\bar{t}$ $(\times 10^2)$", background_dsigma_aa_ttbar, "teal", "--"),
    (r"$\gamma\gamma \to \tau^+\tau^-$", background_dsigma_aa_tautau, "orange", "-"),
    (r"$\gamma\gamma \to \tau^+\tau^- (inel)$", background_dsigma_aa_tautau_inel, "brown", (0, (3, 1, 1, 1))),
    (r"Inclusive $t\bar{t}$", background_dsigma_inclusive_ttbar, "orchid", (0, (3, 1, 1, 1))),
    (r"Single Top", background_dsigma_single_top, "teal", "--"),
    (r"W Production", background_dsigma_w_production, "darkgreen", "-."),
    (r"Z Production", background_dsigma_z_production, "darkred", ":"),
    (r"$WWj$", background_dsigma_wwj_production, "royalblue", "--"),
    (r"$ZZj$ $(\times 10^2)$", background_dsigma_zzj_production, "mediumvioletred", "-."),
    (r"$WZj$", background_dsigma_wzj_production, "slategray", ":")
]

for label, dsigma_dict, color, style in background_styles:
    ht_bins, dsigma = dsigma_dict["ht_total_bins"]
    if sum(dsigma) > 0:
        plt.step(ht_bins, dsigma, where="mid", alpha=0.7,
                 label=label, color=color, linestyle=style, linewidth=3)

# ✅ Axis labels and layout
plt.xlabel(r"$H_T^\mathrm{total} \ \mathrm{[GeV]}$")
plt.ylabel(r"$\frac{d\sigma}{dH_T} \ \mathrm{[pb/GeV]}$")
plt.title(r"Delphes simulation : $e^- p \to e^- W^+ W^- p \to e^- j j \ell \nu_{\ell} p$ : LHeC@1.2 TeV", fontsize=20)
plt.yscale("log")
plt.ylim(1e-5, 1e0)

# ✅ Legend and Save
plt.legend()
plt.grid(True, linestyle="--", alpha=0.6)
plt.tight_layout()
plt.savefig("/home/hamzeh-khanpour/Documents/GitHub/LHeC_Fast_Simulation/Pythia8_Delphes_semi_leptonic_allsignal_bkgs/dsigma_ht_total_allFMsignal_allbkgs.pdf", dpi=600)

plt.show()







# ===================================================
# ✅ Total Scalar Transverse Energy \( H_T \) (Semi-Leptonic)
# ===================================================

plt.figure(figsize=(11, 12))  # New figure for H_T
plt.subplots_adjust(left=0.15, right=0.95, bottom=0.12, top=0.95)

# ✅ Plot Signals
for signal_name, dsigma_data in signal_dsigma.items():
    ht_bins, dsigma = dsigma_data["ht_total_bins"]
    total = sum(dsigma)
    if total > 0:
        dsigma_norm = [x / total for x in dsigma]
        plt.step(ht_bins, dsigma_norm, where="mid", alpha=0.7,
                 label=f"Signal ($W^+W^-$) [{signal_name}]",
                 color=signal_colors.get(signal_name, "black"), linewidth=3)

# ✅ Plot Backgrounds
background_styles = [
    (r"SM ($\gamma\gamma \to W^+W^-$)", background_dsigma_aa_ww, "blue", "-"),
    (r"$\gamma\gamma \to t\bar{t}$ $(\times 10^2)$", background_dsigma_aa_ttbar, "teal", "--"),
    (r"$\gamma\gamma \to \tau^+\tau^-$", background_dsigma_aa_tautau, "orange", "-"),
    (r"$\gamma\gamma \to \tau^+\tau^- (inel)$", background_dsigma_aa_tautau_inel, "brown", (0, (3, 1, 1, 1))),
    (r"Inclusive $t\bar{t}$", background_dsigma_inclusive_ttbar, "orchid", (0, (3, 1, 1, 1))),
    (r"Single Top", background_dsigma_single_top, "teal", "--"),
    (r"W Production", background_dsigma_w_production, "darkgreen", "-."),
    (r"Z Production", background_dsigma_z_production, "darkred", ":"),
    (r"$WWj$", background_dsigma_wwj_production, "royalblue", "--"),
    (r"$ZZj$ $(\times 10^2)$", background_dsigma_zzj_production, "mediumvioletred", "-."),
    (r"$WZj$", background_dsigma_wzj_production, "slategray", ":")
]

for label, dsigma_dict, color, style in background_styles:
    ht_bins, dsigma = dsigma_dict["ht_total_bins"]
    total = sum(dsigma)
    if total > 0:
        dsigma_norm = [x / total for x in dsigma]
        plt.step(ht_bins, dsigma_norm, where="mid", alpha=0.7,
                 label=label, color=color, linestyle=style, linewidth=3)

# ✅ Axis labels and layout
plt.xlabel(r"$H_T^\mathrm{total} \ \mathrm{[GeV]}$")
plt.ylabel(r"Normalized Entries")
plt.title(r"Delphes simulation : $e^- p \to e^- W^+ W^- p \to e^- j j \ell \nu_{\ell} p$ : LHeC@1.2 TeV", fontsize=20)
#plt.yscale("log")
plt.ylim(1e-5, 0.3)

# ✅ Legend and Save
plt.legend()
plt.grid(True, linestyle="--", alpha=0.6)
plt.tight_layout()
plt.savefig("/home/hamzeh-khanpour/Documents/GitHub/LHeC_Fast_Simulation/Pythia8_Delphes_semi_leptonic_allsignal_bkgs/normalized_dsigma_ht_total_allFMsignal_allbkgs.pdf", dpi=600)

plt.show()









# ===================================================
# ✅ Azimuthal Angle Difference Between Jets \( \Delta\phi_{jj} \)
# ===================================================

plt.figure(figsize=(11, 12))  # New figure for Δϕ(jj)
plt.subplots_adjust(left=0.15, right=0.95, bottom=0.12, top=0.95)

# ✅ Plot Signals
for signal_name, dsigma_data in signal_dsigma.items():
    phi_bins, dsigma = dsigma_data["delta_phi_jj_bins"]
    plt.step(phi_bins, dsigma, where="mid", alpha=0.7,
             label=f"Signal ($W^+W^-$) [{signal_name}]",
             color=signal_colors.get(signal_name, "black"), linewidth=3)

# ✅ Background Styles (re-declared to ensure consistency if run standalone)
background_styles = [
    (r"SM ($\gamma\gamma \to W^+W^-$)", background_dsigma_aa_ww, "blue", "-"),
    (r"$\gamma\gamma \to t\bar{t}$ $(\times 10^2)$", background_dsigma_aa_ttbar, "teal", "--"),
    (r"$\gamma\gamma \to \tau^+\tau^-$", background_dsigma_aa_tautau, "orange", "-"),
    (r"$\gamma\gamma \to \tau^+\tau^- (inel)$", background_dsigma_aa_tautau_inel, "brown", (0, (3, 1, 1, 1))),
    (r"Inclusive $t\bar{t}$", background_dsigma_inclusive_ttbar, "orchid", (0, (3, 1, 1, 1))),
    (r"Single Top", background_dsigma_single_top, "teal", "--"),
    (r"W Production", background_dsigma_w_production, "darkgreen", "-."),
    (r"Z Production", background_dsigma_z_production, "darkred", ":"),
    (r"$WWj$", background_dsigma_wwj_production, "royalblue", "--"),
    (r"$ZZj$ $(\times 10^2)$", background_dsigma_zzj_production, "mediumvioletred", "-."),
    (r"$WZj$", background_dsigma_wzj_production, "slategray", ":")
]

# ✅ Plot Backgrounds
for label, dsigma_dict, color, style in background_styles:
    phi_bins, dsigma = dsigma_dict["delta_phi_jj_bins"]
    if sum(dsigma) > 0:
        plt.step(phi_bins, dsigma, where="mid", alpha=0.7,
                 label=label, color=color, linestyle=style, linewidth=3)

# ✅ Axis labels and title
plt.xlabel(r"$\Delta\phi_{jj}$")
plt.ylabel(r"$\frac{d\sigma}{d\Delta\phi_{jj}} \ \mathrm{[pb]}$")
plt.title(r"Delphes simulation : $e^- p \to e^- W^+ W^- p \to e^- j j \ell \nu_{\ell} p$ : LHeC@1.2 TeV", fontsize=20)
plt.yscale("log")
plt.ylim(1e-5, 1e1)

# ✅ Legend, Grid, Save
plt.legend()
plt.grid(True, linestyle="--", alpha=0.6)
plt.tight_layout()
plt.savefig("/home/hamzeh-khanpour/Documents/GitHub/LHeC_Fast_Simulation/Pythia8_Delphes_semi_leptonic_allsignal_bkgs/dsigma_delta_phi_jj_allFMsignal_allbkgs.pdf", dpi=600)

plt.show()






# ===================================================
# ✅ Azimuthal Angle Difference Between Jets \( \Delta\phi_{jj} \)
# ===================================================

plt.figure(figsize=(11, 12))  # New figure for Δϕ(jj)
plt.subplots_adjust(left=0.15, right=0.95, bottom=0.12, top=0.95)

# ✅ Plot Signals
for signal_name, dsigma_data in signal_dsigma.items():
    phi_bins, dsigma = dsigma_data["delta_phi_jj_bins"]
    total = sum(dsigma)
    if total > 0:
        dsigma_norm = [x / total for x in dsigma]
        plt.step(phi_bins, dsigma_norm, where="mid", alpha=0.7,
                 label=f"Signal ($W^+W^-$) [{signal_name}]",
                 color=signal_colors.get(signal_name, "black"), linewidth=3)

# ✅ Background Styles (re-declared to ensure consistency if run standalone)
background_styles = [
    (r"SM ($\gamma\gamma \to W^+W^-$)", background_dsigma_aa_ww, "blue", "-"),
    (r"$\gamma\gamma \to t\bar{t}$ $(\times 10^2)$", background_dsigma_aa_ttbar, "teal", "--"),
    (r"$\gamma\gamma \to \tau^+\tau^-$", background_dsigma_aa_tautau, "orange", "-"),
    (r"$\gamma\gamma \to \tau^+\tau^- (inel)$", background_dsigma_aa_tautau_inel, "brown", (0, (3, 1, 1, 1))),
    (r"Inclusive $t\bar{t}$", background_dsigma_inclusive_ttbar, "orchid", (0, (3, 1, 1, 1))),
    (r"Single Top", background_dsigma_single_top, "teal", "--"),
    (r"W Production", background_dsigma_w_production, "darkgreen", "-."),
    (r"Z Production", background_dsigma_z_production, "darkred", ":"),
    (r"$WWj$", background_dsigma_wwj_production, "royalblue", "--"),
    (r"$ZZj$ $(\times 10^2)$", background_dsigma_zzj_production, "mediumvioletred", "-."),
    (r"$WZj$", background_dsigma_wzj_production, "slategray", ":")
]

# ✅ Plot Backgrounds
for label, dsigma_dict, color, style in background_styles:
    phi_bins, dsigma = dsigma_dict["delta_phi_jj_bins"]
    total = sum(dsigma)
    if total > 0:
        dsigma_norm = [x / total for x in dsigma]
        plt.step(phi_bins, dsigma_norm, where="mid", alpha=0.7,
                 label=label, color=color, linestyle=style, linewidth=3)

# ✅ Axis labels and title
plt.xlabel(r"$\Delta\phi_{jj}$")
plt.ylabel(r"Normalized Entries")
plt.title(r"Delphes simulation : $e^- p \to e^- W^+ W^- p \to e^- j j \ell \nu_{\ell} p$ : LHeC@1.2 TeV", fontsize=20)
#plt.yscale("log")
plt.ylim(1e-5, 0.2)

# ✅ Legend, Grid, Save
plt.legend()
plt.grid(True, linestyle="--", alpha=0.6)
plt.tight_layout()
plt.savefig("/home/hamzeh-khanpour/Documents/GitHub/LHeC_Fast_Simulation/Pythia8_Delphes_semi_leptonic_allsignal_bkgs/normalized_dsigma_delta_phi_jj_allFMsignal_allbkgs.pdf", dpi=600)

plt.show()









# ===================================================
# ✅ Azimuthal Angle Difference Between Leptonic and Hadronic W Bosons \( \Delta\phi_{WW} \)
# ===================================================

plt.figure(figsize=(11, 12))  # New figure for Δϕ(Wlep, Whad)
plt.subplots_adjust(left=0.15, right=0.95, bottom=0.12, top=0.95)

# ✅ Plot Signals
for signal_name, dsigma_data in signal_dsigma.items():
    phi_bins, dsigma = dsigma_data["delta_phi_wl_wh_bins"]
    plt.step(phi_bins, dsigma, where="mid", alpha=0.7,
             label=f"Signal ($W^+W^-$) [{signal_name}]",
             color=signal_colors.get(signal_name, "black"), linewidth=3)

# ✅ Background Styles (ensures standalone use)
background_styles = [
    (r"SM ($\gamma\gamma \to W^+W^-$)", background_dsigma_aa_ww, "blue", "-"),
    (r"$\gamma\gamma \to t\bar{t}$ $(\times 10^2)$", background_dsigma_aa_ttbar, "teal", "--"),
    (r"$\gamma\gamma \to \tau^+\tau^-$", background_dsigma_aa_tautau, "orange", "-"),
    (r"$\gamma\gamma \to \tau^+\tau^- (inel)$", background_dsigma_aa_tautau_inel, "brown", (0, (3, 1, 1, 1))),
    (r"Inclusive $t\bar{t}$", background_dsigma_inclusive_ttbar, "orchid", (0, (3, 1, 1, 1))),
    (r"Single Top", background_dsigma_single_top, "teal", "--"),
    (r"W Production", background_dsigma_w_production, "darkgreen", "-."),
    (r"Z Production", background_dsigma_z_production, "darkred", ":"),
    (r"$WWj$", background_dsigma_wwj_production, "royalblue", "--"),
    (r"$ZZj$ $(\times 10^2)$", background_dsigma_zzj_production, "mediumvioletred", "-."),
    (r"$WZj$", background_dsigma_wzj_production, "slategray", ":")
]

# ✅ Plot Backgrounds
for label, dsigma_dict, color, style in background_styles:
    phi_bins, dsigma = dsigma_dict["delta_phi_wl_wh_bins"]
    if sum(dsigma) > 0:
        plt.step(phi_bins, dsigma, where="mid", alpha=0.7,
                 label=label, color=color, linestyle=style, linewidth=3)

# ✅ Axis Labels and Title
plt.xlabel(r"$\Delta\phi(W^{lep}, W^{had})$")
plt.ylabel(r"$\frac{d\sigma}{d\Delta\phi_{WW}} \ \mathrm{[pb]}$")
plt.title(r"Delphes simulation : $e^- p \to e^- W^+ W^- p \to e^- j j \ell \nu_{\ell} p$ : LHeC@1.2 TeV", fontsize=20)
plt.yscale("log")
plt.ylim(1e-5, 1e0)

# ✅ Legend and Layout
plt.legend()
plt.grid(True, linestyle="--", alpha=0.6)
plt.tight_layout()

# ✅ Save and Show
plt.savefig("/home/hamzeh-khanpour/Documents/GitHub/LHeC_Fast_Simulation/Pythia8_Delphes_semi_leptonic_allsignal_bkgs/dsigma_delta_phi_wl_wh_allFMsignal_allbkgs.pdf", dpi=600)
plt.show()








# ===================================================
# ✅ Azimuthal Angle Difference Between Leptonic and Hadronic W Bosons \( \Delta\phi_{WW} \)
# ===================================================

plt.figure(figsize=(11, 12))  # New figure for Δϕ(Wlep, Whad)
plt.subplots_adjust(left=0.15, right=0.95, bottom=0.12, top=0.95)

# ✅ Plot Signals
for signal_name, dsigma_data in signal_dsigma.items():
    phi_bins, dsigma = dsigma_data["delta_phi_wl_wh_bins"]
    total = sum(dsigma)
    if total > 0:
        dsigma_norm = [x / total for x in dsigma]
        plt.step(phi_bins, dsigma_norm, where="mid", alpha=0.7,
                 label=f"Signal ($W^+W^-$) [{signal_name}]",
                 color=signal_colors.get(signal_name, "black"), linewidth=3)

# ✅ Background Styles (ensures standalone use)
background_styles = [
    (r"SM ($\gamma\gamma \to W^+W^-$)", background_dsigma_aa_ww, "blue", "-"),
    (r"$\gamma\gamma \to t\bar{t}$ $(\times 10^2)$", background_dsigma_aa_ttbar, "teal", "--"),
    (r"$\gamma\gamma \to \tau^+\tau^-$", background_dsigma_aa_tautau, "orange", "-"),
    (r"$\gamma\gamma \to \tau^+\tau^- (inel)$", background_dsigma_aa_tautau_inel, "brown", (0, (3, 1, 1, 1))),
    (r"Inclusive $t\bar{t}$", background_dsigma_inclusive_ttbar, "orchid", (0, (3, 1, 1, 1))),
    (r"Single Top", background_dsigma_single_top, "teal", "--"),
    (r"W Production", background_dsigma_w_production, "darkgreen", "-."),
    (r"Z Production", background_dsigma_z_production, "darkred", ":"),
    (r"$WWj$", background_dsigma_wwj_production, "royalblue", "--"),
    (r"$ZZj$ $(\times 10^2)$", background_dsigma_zzj_production, "mediumvioletred", "-."),
    (r"$WZj$", background_dsigma_wzj_production, "slategray", ":")
]

# ✅ Plot Backgrounds
for label, dsigma_dict, color, style in background_styles:
    phi_bins, dsigma = dsigma_dict["delta_phi_wl_wh_bins"]
    total = sum(dsigma)
    if total > 0:
        dsigma_norm = [x / total for x in dsigma]
        plt.step(phi_bins, dsigma_norm, where="mid", alpha=0.7,
                 label=label, color=color, linestyle=style, linewidth=3)

# ✅ Axis Labels and Title
plt.xlabel(r"$\Delta\phi(W^{lep}, W^{had})$")
plt.ylabel(r"Normalized $\frac{1}{\sigma} \frac{d\sigma}{d\Delta\phi_{WW}}$")
plt.title(r"Delphes simulation : $e^- p \to e^- W^+ W^- p \to e^- j j \ell \nu_{\ell} p$ : LHeC@1.2 TeV", fontsize=20)
#plt.yscale("log")
plt.ylim(1e-5, 1e0)

# ✅ Legend and Layout
plt.legend()
plt.grid(True, linestyle="--", alpha=0.6)
plt.tight_layout()

# ✅ Save and Show
plt.savefig("/home/hamzeh-khanpour/Documents/GitHub/LHeC_Fast_Simulation/Pythia8_Delphes_semi_leptonic_allsignal_bkgs/normalized_dsigma_delta_phi_wl_wh_allFMsignal_allbkgs.pdf", dpi=600)
plt.show()










# ===================================================
# ✅ Pseudorapidity Difference Between Leptonic and Hadronic W Bosons \( \Delta\eta_{WW} \)
# ===================================================

plt.figure(figsize=(11, 12))  # New figure for Δη(Wlep, Whad)
plt.subplots_adjust(left=0.15, right=0.95, bottom=0.12, top=0.95)

# ✅ Plot Signals
for signal_name, dsigma_data in signal_dsigma.items():
    eta_bins, dsigma = dsigma_data["delta_eta_wl_wh_bins"]
    plt.step(eta_bins, dsigma, where="mid", alpha=0.7,
             label=f"Signal ($W^+W^-$) [{signal_name}]",
             color=signal_colors.get(signal_name, "black"), linewidth=3)

# ✅ Background Styles
background_styles = [
    (r"SM ($\gamma\gamma \to W^+W^-$)", background_dsigma_aa_ww, "blue", "-"),
    (r"$\gamma\gamma \to t\bar{t}$ $(\times 10^2)$", background_dsigma_aa_ttbar, "teal", "--"),
    (r"$\gamma\gamma \to \tau^+\tau^-$", background_dsigma_aa_tautau, "orange", "-"),
    (r"$\gamma\gamma \to \tau^+\tau^- (inel)$", background_dsigma_aa_tautau_inel, "brown", (0, (3, 1, 1, 1))),
    (r"Inclusive $t\bar{t}$", background_dsigma_inclusive_ttbar, "orchid", (0, (3, 1, 1, 1))),
    (r"Single Top", background_dsigma_single_top, "teal", "--"),
    (r"W Production", background_dsigma_w_production, "darkgreen", "-."),
    (r"Z Production", background_dsigma_z_production, "darkred", ":"),
    (r"$WWj$", background_dsigma_wwj_production, "royalblue", "--"),
    (r"$ZZj$ $(\times 10^2)$", background_dsigma_zzj_production, "mediumvioletred", "-."),
    (r"$WZj$", background_dsigma_wzj_production, "slategray", ":")
]

# ✅ Plot Backgrounds
for label, dsigma_dict, color, style in background_styles:
    eta_bins, dsigma = dsigma_dict["delta_eta_wl_wh_bins"]
    if sum(dsigma) > 0:
        plt.step(eta_bins, dsigma, where="mid", alpha=0.7,
                 label=label, color=color, linestyle=style, linewidth=3)

# ✅ Axis Labels and Title
plt.xlabel(r"$\Delta\eta(W^{lep}, W^{had})$")
plt.ylabel(r"$\frac{d\sigma}{d\Delta\eta_{WW}} \ \mathrm{[pb]}$")
plt.title(r"Delphes simulation : $e^- p \to e^- W^+ W^- p \to e^- j j \ell \nu_{\ell} p$ : LHeC@1.2 TeV", fontsize=20)
plt.yscale("log")
plt.ylim(1e-5, 1e1)

# ✅ Legend, Grid, Save
plt.legend()
plt.grid(True, linestyle="--", alpha=0.6)
plt.tight_layout()
plt.savefig("/home/hamzeh-khanpour/Documents/GitHub/LHeC_Fast_Simulation/Pythia8_Delphes_semi_leptonic_allsignal_bkgs/dsigma_delta_eta_wl_wh_allFMsignal_allbkgs.pdf", dpi=600)
plt.show()










# ===================================================
# ✅ Pseudorapidity Difference Between Leptonic and Hadronic W Bosons \( \Delta\eta_{WW} \)
# ===================================================

plt.figure(figsize=(11, 12))  # New figure for Δη(Wlep, Whad)
plt.subplots_adjust(left=0.15, right=0.95, bottom=0.12, top=0.95)

# ✅ Plot Signals
for signal_name, dsigma_data in signal_dsigma.items():
    eta_bins, dsigma = dsigma_data["delta_eta_wl_wh_bins"]
    total = sum(dsigma)
    if total > 0:
        dsigma_norm = [x / total for x in dsigma]
        plt.step(eta_bins, dsigma_norm, where="mid", alpha=0.7,
                 label=f"Signal ($W^+W^-$) [{signal_name}]",
                 color=signal_colors.get(signal_name, "black"), linewidth=3)

# ✅ Background Styles
background_styles = [
    (r"SM ($\gamma\gamma \to W^+W^-$)", background_dsigma_aa_ww, "blue", "-"),
    (r"$\gamma\gamma \to t\bar{t}$ $(\times 10^2)$", background_dsigma_aa_ttbar, "teal", "--"),
    (r"$\gamma\gamma \to \tau^+\tau^-$", background_dsigma_aa_tautau, "orange", "-"),
    (r"$\gamma\gamma \to \tau^+\tau^- (inel)$", background_dsigma_aa_tautau_inel, "brown", (0, (3, 1, 1, 1))),
    (r"Inclusive $t\bar{t}$", background_dsigma_inclusive_ttbar, "orchid", (0, (3, 1, 1, 1))),
    (r"Single Top", background_dsigma_single_top, "teal", "--"),
    (r"W Production", background_dsigma_w_production, "darkgreen", "-."),
    (r"Z Production", background_dsigma_z_production, "darkred", ":"),
    (r"$WWj$", background_dsigma_wwj_production, "royalblue", "--"),
    (r"$ZZj$ $(\times 10^2)$", background_dsigma_zzj_production, "mediumvioletred", "-."),
    (r"$WZj$", background_dsigma_wzj_production, "slategray", ":")
]

# ✅ Plot Backgrounds
for label, dsigma_dict, color, style in background_styles:
    eta_bins, dsigma = dsigma_dict["delta_eta_wl_wh_bins"]
    total = sum(dsigma)
    if total > 0:
        dsigma_norm = [x / total for x in dsigma]
        plt.step(eta_bins, dsigma_norm, where="mid", alpha=0.7,
                 label=label, color=color, linestyle=style, linewidth=3)

# ✅ Axis Labels and Title
plt.xlabel(r"$\Delta\eta(W^{lep}, W^{had})$")
plt.ylabel(r"Normalized $\frac{1}{\sigma} \frac{d\sigma}{d\Delta\eta_{WW}}$")
plt.title(r"Delphes simulation : $e^- p \to e^- W^+ W^- p \to e^- j j \ell \nu_{\ell} p$ : LHeC@1.2 TeV", fontsize=20)
#plt.yscale("log")
plt.ylim(1e-5, 0.2)

# ✅ Legend, Grid, Save
plt.legend()
plt.grid(True, linestyle="--", alpha=0.6)
plt.tight_layout()
plt.savefig("/home/hamzeh-khanpour/Documents/GitHub/LHeC_Fast_Simulation/Pythia8_Delphes_semi_leptonic_allsignal_bkgs/normalized_dsigma_delta_eta_wl_wh_allFMsignal_allbkgs.pdf", dpi=600)
plt.show()









# ===================================================
# ✅ Dijet Invariant Mass \( m_{jj} \) Differential Cross-Section
# ===================================================

plt.figure(figsize=(11, 12))  # New figure for m_jj
plt.subplots_adjust(left=0.15, right=0.95, bottom=0.12, top=0.95)

# ✅ Plot Signals
for signal_name, dsigma_data in signal_dsigma.items():
    mjj_bins, dsigma = dsigma_data["m_jj_bins"]
    plt.step(mjj_bins, dsigma, where="mid", alpha=0.7,
             label=f"Signal ($W^+W^-$) [{signal_name}]",
             color=signal_colors.get(signal_name, "black"), linewidth=3)

# ✅ Background styles (reuse consistent style block if not global)
background_styles = [
    (r"SM ($\gamma\gamma \to W^+W^-$)", background_dsigma_aa_ww, "blue", "-"),
    (r"$\gamma\gamma \to t\bar{t}$ $(\times 10^2)$", background_dsigma_aa_ttbar, "teal", "--"),
    (r"$\gamma\gamma \to \tau^+\tau^-$", background_dsigma_aa_tautau, "orange", "-"),
    (r"$\gamma\gamma \to \tau^+\tau^- (inel)$", background_dsigma_aa_tautau_inel, "brown", (0, (3, 1, 1, 1))),
    (r"Inclusive $t\bar{t}$", background_dsigma_inclusive_ttbar, "orchid", (0, (3, 1, 1, 1))),
    (r"Single Top", background_dsigma_single_top, "teal", "--"),
    (r"W Production", background_dsigma_w_production, "darkgreen", "-."),
    (r"Z Production", background_dsigma_z_production, "darkred", ":"),
    (r"$WWj$", background_dsigma_wwj_production, "royalblue", "--"),
    (r"$ZZj$ $(\times 10^2)$", background_dsigma_zzj_production, "mediumvioletred", "-."),
    (r"$WZj$", background_dsigma_wzj_production, "slategray", ":")
]

# ✅ Plot Backgrounds
for label, dsigma_dict, color, style in background_styles:
    mjj_bins, dsigma = dsigma_dict["m_jj_bins"]
    if sum(dsigma) > 0:
        plt.step(mjj_bins, dsigma, where="mid", alpha=0.7,
                 label=label, color=color, linestyle=style, linewidth=3)

# ✅ Axis labels and title
plt.xlabel(r"$m_{jj} \ \mathrm{[GeV]}$")
plt.ylabel(r"$\frac{d\sigma}{dm_{jj}} \ \mathrm{[pb/GeV]}$")
plt.title(r"Delphes simulation : $e^- p \to e^- W^+ W^- p \to e^- j j \ell \nu_{\ell} p$ : LHeC@1.2 TeV", fontsize=20)
plt.yscale("log")
plt.ylim(1e-5, 1e0)

# ✅ Legend and Save
plt.legend()
plt.grid(True, linestyle="--", alpha=0.6)
plt.tight_layout()
plt.savefig("/home/hamzeh-khanpour/Documents/GitHub/LHeC_Fast_Simulation/Pythia8_Delphes_semi_leptonic_allsignal_bkgs/dsigma_m_jj_allFMsignal_allbkgs.pdf", dpi=600)

plt.show()








# ===================================================
# ✅ Dijet Invariant Mass \( m_{jj} \) Differential Cross-Section (Normalized)
# ===================================================

plt.figure(figsize=(11, 12))  # New figure for m_jj
plt.subplots_adjust(left=0.15, right=0.95, bottom=0.12, top=0.95)

# ✅ Plot Signals
for signal_name, dsigma_data in signal_dsigma.items():
    mjj_bins, dsigma = dsigma_data["m_jj_bins"]
    total = sum(dsigma)
    if total > 0:
        dsigma_norm = [x / total for x in dsigma]
        plt.step(mjj_bins, dsigma_norm, where="mid", alpha=0.7,
                 label=f"Signal ($W^+W^-$) [{signal_name}]",
                 color=signal_colors.get(signal_name, "black"), linewidth=3)

# ✅ Background styles (reuse consistent style block if not global)
background_styles = [
    (r"SM ($\gamma\gamma \to W^+W^-$)", background_dsigma_aa_ww, "blue", "-"),
    (r"$\gamma\gamma \to t\bar{t}$ $(\times 10^2)$", background_dsigma_aa_ttbar, "teal", "--"),
    (r"$\gamma\gamma \to \tau^+\tau^-$", background_dsigma_aa_tautau, "orange", "-"),
    (r"$\gamma\gamma \to \tau^+\tau^- (inel)$", background_dsigma_aa_tautau_inel, "brown", (0, (3, 1, 1, 1))),
    (r"Inclusive $t\bar{t}$", background_dsigma_inclusive_ttbar, "orchid", (0, (3, 1, 1, 1))),
    (r"Single Top", background_dsigma_single_top, "teal", "--"),
    (r"W Production", background_dsigma_w_production, "darkgreen", "-."),
    (r"Z Production", background_dsigma_z_production, "darkred", ":"),
    (r"$WWj$", background_dsigma_wwj_production, "royalblue", "--"),
    (r"$ZZj$ $(\times 10^2)$", background_dsigma_zzj_production, "mediumvioletred", "-."),
    (r"$WZj$", background_dsigma_wzj_production, "slategray", ":")
]

# ✅ Plot Backgrounds
for label, dsigma_dict, color, style in background_styles:
    mjj_bins, dsigma = dsigma_dict["m_jj_bins"]
    total = sum(dsigma)
    if total > 0:
        dsigma_norm = [x / total for x in dsigma]
        plt.step(mjj_bins, dsigma_norm, where="mid", alpha=0.7,
                 label=label, color=color, linestyle=style, linewidth=3)

# ✅ Axis labels and title
plt.xlabel(r"$m_{jj} \ \mathrm{[GeV]}$")
plt.ylabel(r"Normalized Entries")
plt.title(r"Delphes simulation : $e^- p \to e^- W^+ W^- p \to e^- j j \ell \nu_{\ell} p$ : LHeC@1.2 TeV", fontsize=20)
#plt.yscale("log")
plt.ylim(1e-5, 0.6)

# ✅ Legend and Save
plt.legend()
plt.grid(True, linestyle="--", alpha=0.6)
plt.tight_layout()
plt.savefig("/home/hamzeh-khanpour/Documents/GitHub/LHeC_Fast_Simulation/Pythia8_Delphes_semi_leptonic_allsignal_bkgs/normalized_dsigma_m_jj_allFMsignal_allbkgs.pdf", dpi=600)

plt.show()










# ===================================================
# ✅ Invariant Mass \( m(\ell\nu jj) \) Differential Cross-Section
# ===================================================

plt.figure(figsize=(11, 12))  # New figure for m_lvjj
plt.subplots_adjust(left=0.15, right=0.95, bottom=0.12, top=0.95)

# ✅ Plot Signals
for signal_name, dsigma_data in signal_dsigma.items():
    mlvjj_bins, dsigma = dsigma_data["m_lvjj_bins"]
    plt.step(mlvjj_bins, dsigma, where="mid", alpha=0.7,
             label=f"Signal ($W^+W^-$) [{signal_name}]",
             color=signal_colors.get(signal_name, "black"), linewidth=3)

# ✅ Backgrounds (same consistent style)
background_styles = [
    (r"SM ($\gamma\gamma \to W^+W^-$)", background_dsigma_aa_ww, "blue", "-"),
    (r"$\gamma\gamma \to t\bar{t}$ $(\times 10^2)$", background_dsigma_aa_ttbar, "teal", "--"),
    (r"$\gamma\gamma \to \tau^+\tau^-$", background_dsigma_aa_tautau, "orange", "-"),
    (r"$\gamma\gamma \to \tau^+\tau^- (inel)$", background_dsigma_aa_tautau_inel, "brown", (0, (3, 1, 1, 1))),
    (r"Inclusive $t\bar{t}$", background_dsigma_inclusive_ttbar, "orchid", (0, (3, 1, 1, 1))),
    (r"Single Top", background_dsigma_single_top, "teal", "--"),
    (r"W Production", background_dsigma_w_production, "darkgreen", "-."),
    (r"Z Production", background_dsigma_z_production, "darkred", ":"),
    (r"$WWj$", background_dsigma_wwj_production, "royalblue", "--"),
    (r"$ZZj$ $(\times 10^2)$", background_dsigma_zzj_production, "mediumvioletred", "-."),
    (r"$WZj$", background_dsigma_wzj_production, "slategray", ":")
]

# ✅ Plot Backgrounds
for label, dsigma_dict, color, style in background_styles:
    mlvjj_bins, dsigma = dsigma_dict["m_lvjj_bins"]
    if sum(dsigma) > 0:
        plt.step(mlvjj_bins, dsigma, where="mid", alpha=0.7,
                 label=label, color=color, linestyle=style, linewidth=3)

# ✅ Axis labels and title
plt.xlabel(r"$m(\ell \nu jj) \ \mathrm{[GeV]}$")
plt.ylabel(r"$\frac{d\sigma}{dm(\ell \nu jj)} \ \mathrm{[pb/GeV]}$")
plt.title(r"Delphes simulation : $e^- p \to e^- W^+ W^- p \to e^- j j \ell \nu_{\ell} p$ : LHeC@1.2 TeV", fontsize=20)
plt.yscale("log")
plt.ylim(1e-5, 1e0)

# ✅ Legend and Save
plt.legend()
plt.grid(True, linestyle="--", alpha=0.6)
plt.tight_layout()
plt.savefig("/home/hamzeh-khanpour/Documents/GitHub/LHeC_Fast_Simulation/Pythia8_Delphes_semi_leptonic_allsignal_bkgs/dsigma_m_lvjj_allFMsignal_allbkgs.pdf", dpi=600)

plt.show()





# ===================================================
# ✅ Invariant Mass \( m(\ell\nu jj) \) Differential Cross-Section (Normalized)
# ===================================================

plt.figure(figsize=(11, 12))  # New figure for m_lvjj
plt.subplots_adjust(left=0.15, right=0.95, bottom=0.12, top=0.95)

# ✅ Plot Signals (normalized)
for signal_name, dsigma_data in signal_dsigma.items():
    mlvjj_bins, dsigma = dsigma_data["m_lvjj_bins"]
    total = sum(dsigma)
    if total > 0:
        dsigma_norm = [x / total for x in dsigma]
        plt.step(mlvjj_bins, dsigma_norm, where="mid", alpha=0.7,
                 label=f"Signal ($W^+W^-$) [{signal_name}]",
                 color=signal_colors.get(signal_name, "black"), linewidth=3)

# ✅ Backgrounds (same consistent style, normalized)
background_styles = [
    (r"SM ($\gamma\gamma \to W^+W^-$)", background_dsigma_aa_ww, "blue", "-"),
    (r"$\gamma\gamma \to t\bar{t}$ $(\times 10^2)$", background_dsigma_aa_ttbar, "teal", "--"),
    (r"$\gamma\gamma \to \tau^+\tau^-$", background_dsigma_aa_tautau, "orange", "-"),
    (r"$\gamma\gamma \to \tau^+\tau^- (inel)$", background_dsigma_aa_tautau_inel, "brown", (0, (3, 1, 1, 1))),
    (r"Inclusive $t\bar{t}$", background_dsigma_inclusive_ttbar, "orchid", (0, (3, 1, 1, 1))),
    (r"Single Top", background_dsigma_single_top, "teal", "--"),
    (r"W Production", background_dsigma_w_production, "darkgreen", "-."),
    (r"Z Production", background_dsigma_z_production, "darkred", ":"),
    (r"$WWj$", background_dsigma_wwj_production, "royalblue", "--"),
    (r"$ZZj$ $(\times 10^2)$", background_dsigma_zzj_production, "mediumvioletred", "-."),
    (r"$WZj$", background_dsigma_wzj_production, "slategray", ":")
]

# ✅ Plot Backgrounds (normalized)
for label, dsigma_dict, color, style in background_styles:
    mlvjj_bins, dsigma = dsigma_dict["m_lvjj_bins"]
    total = sum(dsigma)
    if total > 0:
        dsigma_norm = [x / total for x in dsigma]
        plt.step(mlvjj_bins, dsigma_norm, where="mid", alpha=0.7,
                 label=label, color=color, linestyle=style, linewidth=3)

# ✅ Axis labels and title
plt.xlabel(r"$m(\ell \nu jj) \ \mathrm{[GeV]}$")
plt.ylabel(r"Normalized Entries")
plt.title(r"Delphes simulation : $e^- p \to e^- W^+ W^- p \to e^- j j \ell \nu_{\ell} p$ : LHeC@1.2 TeV", fontsize=20)
#plt.yscale("log")
plt.ylim(1e-5, 0.2)

# ✅ Legend and Save
plt.legend()
plt.grid(True, linestyle="--", alpha=0.6)
plt.tight_layout()
plt.savefig("/home/hamzeh-khanpour/Documents/GitHub/LHeC_Fast_Simulation/Pythia8_Delphes_semi_leptonic_allsignal_bkgs/normalized_dsigma_m_lvjj_allFMsignal_allbkgs.pdf", dpi=600)

plt.show()








#=========================================================================
#=========================================================================

import ROOT

# ✅ Set the output ROOT file path
output_file = ROOT.TFile("/home/hamzeh-khanpour/Delphes-3.5.0/output_histograms_FM1.root", "RECREATE")

try:
    # ✅ Create signal directories with sanitized names
    signal_dirs = {}
    for signal_name in signal_files.keys():
        signal_dirs[signal_name] = output_file.mkdir(f"signal_{signal_name}")


    background_histogram_sets = {
        "aa_ww": {
            "hist_lepton_pt": hist_lepton_pt_aa_ww,
            "hist_leading_jet_pt": hist_leading_jet_pt_aa_ww,
            "hist_lepton_eta": hist_lepton_eta_aa_ww,
            "hist_delta_r": hist_delta_r_aa_ww,
            "hist_missing_et": hist_missing_et_aa_ww,
            "hist_subleading_jet_eta": hist_subleading_jet_eta_aa_ww,
            "hist_leading_jet_eta": hist_leading_jet_eta_aa_ww,
            "hist_jet_centrality": hist_jet_centrality_aa_ww,
            "hist_delta_eta_jj": hist_delta_eta_jj_aa_ww,
            "hist_m_w_hadronic": hist_m_w_hadronic_aa_ww,
            "hist_m_w_leptonic": hist_m_w_leptonic_aa_ww,
            # ✅ Newly added histograms:
            "hist_pt_w_leptonic": hist_pt_w_leptonic_aa_ww,
            "hist_pt_w_hadronic": hist_pt_w_hadronic_aa_ww,
            "hist_delta_phi_lep_met": hist_delta_phi_lep_met_aa_ww,
            "hist_mt_w_leptonic": hist_mt_w_leptonic_aa_ww,
            "hist_ht_total": hist_ht_total_aa_ww,
            "hist_delta_phi_jj": hist_delta_phi_jj_aa_ww,
            "hist_delta_phi_wl_wh": hist_delta_phi_wl_wh_aa_ww,
            "hist_delta_eta_wl_wh": hist_delta_eta_wl_wh_aa_ww,
            "hist_m_jj": hist_m_jj_aa_ww,
            "hist_m_lvjj": hist_m_lvjj_aa_ww
        },
        "aa_ttbar": {
            "hist_lepton_pt": hist_lepton_pt_aa_ttbar,
            "hist_leading_jet_pt": hist_leading_jet_pt_aa_ttbar,
            "hist_lepton_eta": hist_lepton_eta_aa_ttbar,
            "hist_delta_r": hist_delta_r_aa_ttbar,
            "hist_missing_et": hist_missing_et_aa_ttbar,
            "hist_subleading_jet_eta": hist_subleading_jet_eta_aa_ttbar,
            "hist_leading_jet_eta": hist_leading_jet_eta_aa_ttbar,
            "hist_jet_centrality": hist_jet_centrality_aa_ttbar,
            "hist_delta_eta_jj": hist_delta_eta_jj_aa_ttbar,
            "hist_m_w_hadronic": hist_m_w_hadronic_aa_ttbar,
            "hist_m_w_leptonic": hist_m_w_leptonic_aa_ttbar,
            # ✅ Newly added observables:
            "hist_pt_w_leptonic": hist_pt_w_leptonic_aa_ttbar,
            "hist_pt_w_hadronic": hist_pt_w_hadronic_aa_ttbar,
            "hist_delta_phi_lep_met": hist_delta_phi_lep_met_aa_ttbar,
            "hist_mt_w_leptonic": hist_mt_w_leptonic_aa_ttbar,
            "hist_ht_total": hist_ht_total_aa_ttbar,
            "hist_delta_phi_jj": hist_delta_phi_jj_aa_ttbar,
            "hist_delta_phi_wl_wh": hist_delta_phi_wl_wh_aa_ttbar,
            "hist_delta_eta_wl_wh": hist_delta_eta_wl_wh_aa_ttbar,
            "hist_m_jj": hist_m_jj_aa_ttbar,
            "hist_m_lvjj": hist_m_lvjj_aa_ttbar
        },
        "aa_tautau": {
            "hist_lepton_pt": hist_lepton_pt_aa_tautau,
            "hist_leading_jet_pt": hist_leading_jet_pt_aa_tautau,
            "hist_lepton_eta": hist_lepton_eta_aa_tautau,
            "hist_delta_r": hist_delta_r_aa_tautau,
            "hist_missing_et": hist_missing_et_aa_tautau,
            "hist_subleading_jet_eta": hist_subleading_jet_eta_aa_tautau,
            "hist_leading_jet_eta": hist_leading_jet_eta_aa_tautau,
            "hist_jet_centrality": hist_jet_centrality_aa_tautau,
            "hist_delta_eta_jj": hist_delta_eta_jj_aa_tautau,
            "hist_m_w_hadronic": hist_m_w_hadronic_aa_tautau,
            "hist_m_w_leptonic": hist_m_w_leptonic_aa_tautau,
            # ✅ New variables
            "hist_pt_w_leptonic": hist_pt_w_leptonic_aa_tautau,
            "hist_pt_w_hadronic": hist_pt_w_hadronic_aa_tautau,
            "hist_delta_phi_lep_met": hist_delta_phi_lep_met_aa_tautau,
            "hist_mt_w_leptonic": hist_mt_w_leptonic_aa_tautau,
            "hist_ht_total": hist_ht_total_aa_tautau,
            "hist_delta_phi_jj": hist_delta_phi_jj_aa_tautau,
            "hist_delta_phi_wl_wh": hist_delta_phi_wl_wh_aa_tautau,
            "hist_delta_eta_wl_wh": hist_delta_eta_wl_wh_aa_tautau,
            "hist_m_jj": hist_m_jj_aa_tautau,
            "hist_m_lvjj": hist_m_lvjj_aa_tautau
        },
        "aa_tautau_inel": {
            "hist_lepton_pt": hist_lepton_pt_aa_tautau_inel,
            "hist_leading_jet_pt": hist_leading_jet_pt_aa_tautau_inel,
            "hist_lepton_eta": hist_lepton_eta_aa_tautau_inel,
            "hist_delta_r": hist_delta_r_aa_tautau_inel,
            "hist_missing_et": hist_missing_et_aa_tautau_inel,
            "hist_subleading_jet_eta": hist_subleading_jet_eta_aa_tautau_inel,
            "hist_leading_jet_eta": hist_leading_jet_eta_aa_tautau_inel,
            "hist_jet_centrality": hist_jet_centrality_aa_tautau_inel,
            "hist_delta_eta_jj": hist_delta_eta_jj_aa_tautau_inel,
            "hist_m_w_hadronic": hist_m_w_hadronic_aa_tautau_inel,
            "hist_m_w_leptonic": hist_m_w_leptonic_aa_tautau_inel,
            # ✅ New observables
            "hist_pt_w_leptonic": hist_pt_w_leptonic_aa_tautau_inel,
            "hist_pt_w_hadronic": hist_pt_w_hadronic_aa_tautau_inel,
            "hist_delta_phi_lep_met": hist_delta_phi_lep_met_aa_tautau_inel,
            "hist_mt_w_leptonic": hist_mt_w_leptonic_aa_tautau_inel,
            "hist_ht_total": hist_ht_total_aa_tautau_inel,
            "hist_delta_phi_jj": hist_delta_phi_jj_aa_tautau_inel,
            "hist_delta_phi_wl_wh": hist_delta_phi_wl_wh_aa_tautau_inel,
            "hist_delta_eta_wl_wh": hist_delta_eta_wl_wh_aa_tautau_inel,
            "hist_m_jj": hist_m_jj_aa_tautau_inel,
            "hist_m_lvjj": hist_m_lvjj_aa_tautau_inel
        },
        "inclusive_ttbar": {
            "hist_lepton_pt": hist_lepton_pt_inclusive_ttbar,
            "hist_leading_jet_pt": hist_leading_jet_pt_inclusive_ttbar,
            "hist_lepton_eta": hist_lepton_eta_inclusive_ttbar,
            "hist_delta_r": hist_delta_r_inclusive_ttbar,
            "hist_missing_et": hist_missing_et_inclusive_ttbar,
            "hist_subleading_jet_eta": hist_subleading_jet_eta_inclusive_ttbar,
            "hist_leading_jet_eta": hist_leading_jet_eta_inclusive_ttbar,
            "hist_jet_centrality": hist_jet_centrality_inclusive_ttbar,
            "hist_delta_eta_jj": hist_delta_eta_jj_inclusive_ttbar,
            "hist_m_w_hadronic": hist_m_w_hadronic_inclusive_ttbar,
            "hist_m_w_leptonic": hist_m_w_leptonic_inclusive_ttbar,
            # ✅ New observables
            "hist_pt_w_leptonic": hist_pt_w_leptonic_inclusive_ttbar,
            "hist_pt_w_hadronic": hist_pt_w_hadronic_inclusive_ttbar,
            "hist_delta_phi_lep_met": hist_delta_phi_lep_met_inclusive_ttbar,
            "hist_mt_w_leptonic": hist_mt_w_leptonic_inclusive_ttbar,
            "hist_ht_total": hist_ht_total_inclusive_ttbar,
            "hist_delta_phi_jj": hist_delta_phi_jj_inclusive_ttbar,
            "hist_delta_phi_wl_wh": hist_delta_phi_wl_wh_inclusive_ttbar,
            "hist_delta_eta_wl_wh": hist_delta_eta_wl_wh_inclusive_ttbar,
            "hist_m_jj": hist_m_jj_inclusive_ttbar,
            "hist_m_lvjj": hist_m_lvjj_inclusive_ttbar
        },
        "single_top": {
            "hist_lepton_pt": hist_lepton_pt_single_top,
            "hist_leading_jet_pt": hist_leading_jet_pt_single_top,
            "hist_lepton_eta": hist_lepton_eta_single_top,
            "hist_delta_r": hist_delta_r_single_top,
            "hist_missing_et": hist_missing_et_single_top,
            "hist_subleading_jet_eta": hist_subleading_jet_eta_single_top,
            "hist_leading_jet_eta": hist_leading_jet_eta_single_top,
            "hist_jet_centrality": hist_jet_centrality_single_top,
            "hist_delta_eta_jj": hist_delta_eta_jj_single_top,
            "hist_m_w_hadronic": hist_m_w_hadronic_single_top,
            "hist_m_w_leptonic": hist_m_w_leptonic_single_top,
            # ✅ New observables
            "hist_pt_w_leptonic": hist_pt_w_leptonic_single_top,
            "hist_pt_w_hadronic": hist_pt_w_hadronic_single_top,
            "hist_delta_phi_lep_met": hist_delta_phi_lep_met_single_top,
            "hist_mt_w_leptonic": hist_mt_w_leptonic_single_top,
            "hist_ht_total": hist_ht_total_single_top,
            "hist_delta_phi_jj": hist_delta_phi_jj_single_top,
            "hist_delta_phi_wl_wh": hist_delta_phi_wl_wh_single_top,
            "hist_delta_eta_wl_wh": hist_delta_eta_wl_wh_single_top,
            "hist_m_jj": hist_m_jj_single_top,
            "hist_m_lvjj": hist_m_lvjj_single_top
        },
        "w_production": {
            "hist_lepton_pt": hist_lepton_pt_w_production,
            "hist_leading_jet_pt": hist_leading_jet_pt_w_production,
            "hist_lepton_eta": hist_lepton_eta_w_production,
            "hist_delta_r": hist_delta_r_w_production,
            "hist_missing_et": hist_missing_et_w_production,
            "hist_subleading_jet_eta": hist_subleading_jet_eta_w_production,
            "hist_leading_jet_eta": hist_leading_jet_eta_w_production,
            "hist_jet_centrality": hist_jet_centrality_w_production,
            "hist_delta_eta_jj": hist_delta_eta_jj_w_production,
            "hist_m_w_hadronic": hist_m_w_hadronic_w_production,
            "hist_m_w_leptonic": hist_m_w_leptonic_w_production,
            # ✅ New observables
            "hist_pt_w_leptonic": hist_pt_w_leptonic_w_production,
            "hist_pt_w_hadronic": hist_pt_w_hadronic_w_production,
            "hist_delta_phi_lep_met": hist_delta_phi_lep_met_w_production,
            "hist_mt_w_leptonic": hist_mt_w_leptonic_w_production,
            "hist_ht_total": hist_ht_total_w_production,
            "hist_delta_phi_jj": hist_delta_phi_jj_w_production,
            "hist_delta_phi_wl_wh": hist_delta_phi_wl_wh_w_production,
            "hist_delta_eta_wl_wh": hist_delta_eta_wl_wh_w_production,
            "hist_m_jj": hist_m_jj_w_production,
            "hist_m_lvjj": hist_m_lvjj_w_production
        },
        "z_production": {
            "hist_lepton_pt": hist_lepton_pt_z_production,
            "hist_leading_jet_pt": hist_leading_jet_pt_z_production,
            "hist_lepton_eta": hist_lepton_eta_z_production,
            "hist_delta_r": hist_delta_r_z_production,
            "hist_missing_et": hist_missing_et_z_production,
            "hist_subleading_jet_eta": hist_subleading_jet_eta_z_production,
            "hist_leading_jet_eta": hist_leading_jet_eta_z_production,
            "hist_jet_centrality": hist_jet_centrality_z_production,
            "hist_delta_eta_jj": hist_delta_eta_jj_z_production,
            "hist_m_w_hadronic": hist_m_w_hadronic_z_production,
            "hist_m_w_leptonic": hist_m_w_leptonic_z_production,
            # ✅ New observables
            "hist_pt_w_leptonic": hist_pt_w_leptonic_z_production,
            "hist_pt_w_hadronic": hist_pt_w_hadronic_z_production,
            "hist_delta_phi_lep_met": hist_delta_phi_lep_met_z_production,
            "hist_mt_w_leptonic": hist_mt_w_leptonic_z_production,
            "hist_ht_total": hist_ht_total_z_production,
            "hist_delta_phi_jj": hist_delta_phi_jj_z_production,
            "hist_delta_phi_wl_wh": hist_delta_phi_wl_wh_z_production,
            "hist_delta_eta_wl_wh": hist_delta_eta_wl_wh_z_production,
            "hist_m_jj": hist_m_jj_z_production,
            "hist_m_lvjj": hist_m_lvjj_z_production
        },
        "wwj": {
            "hist_lepton_pt": hist_lepton_pt_wwj_production,
            "hist_leading_jet_pt": hist_leading_jet_pt_wwj_production,
            "hist_lepton_eta": hist_lepton_eta_wwj_production,
            "hist_delta_r": hist_delta_r_wwj_production,
            "hist_missing_et": hist_missing_et_wwj_production,
            "hist_subleading_jet_eta": hist_subleading_jet_eta_wwj_production,
            "hist_leading_jet_eta": hist_leading_jet_eta_wwj_production,
            "hist_jet_centrality": hist_jet_centrality_wwj_production,
            "hist_delta_eta_jj": hist_delta_eta_jj_wwj_production,
            "hist_m_w_hadronic": hist_m_w_hadronic_wwj_production,
            "hist_m_w_leptonic": hist_m_w_leptonic_wwj_production,
            # ✅ New observables
            "hist_pt_w_leptonic": hist_pt_w_leptonic_wwj_production,
            "hist_pt_w_hadronic": hist_pt_w_hadronic_wwj_production,
            "hist_delta_phi_lep_met": hist_delta_phi_lep_met_wwj_production,
            "hist_mt_w_leptonic": hist_mt_w_leptonic_wwj_production,
            "hist_ht_total": hist_ht_total_wwj_production,
            "hist_delta_phi_jj": hist_delta_phi_jj_wwj_production,
            "hist_delta_phi_wl_wh": hist_delta_phi_wl_wh_wwj_production,
            "hist_delta_eta_wl_wh": hist_delta_eta_wl_wh_wwj_production,
            "hist_m_jj": hist_m_jj_wwj_production,
            "hist_m_lvjj": hist_m_lvjj_wwj_production
        },
        "zzj": {
            "hist_lepton_pt": hist_lepton_pt_zzj_production,
            "hist_leading_jet_pt": hist_leading_jet_pt_zzj_production,
            "hist_lepton_eta": hist_lepton_eta_zzj_production,
            "hist_delta_r": hist_delta_r_zzj_production,
            "hist_missing_et": hist_missing_et_zzj_production,
            "hist_subleading_jet_eta": hist_subleading_jet_eta_zzj_production,
            "hist_leading_jet_eta": hist_leading_jet_eta_zzj_production,
            "hist_jet_centrality": hist_jet_centrality_zzj_production,
            "hist_delta_eta_jj": hist_delta_eta_jj_zzj_production,
            "hist_m_w_hadronic": hist_m_w_hadronic_zzj_production,
            "hist_m_w_leptonic": hist_m_w_leptonic_zzj_production,
            # ✅ New variables
            "hist_pt_w_leptonic": hist_pt_w_leptonic_zzj_production,
            "hist_pt_w_hadronic": hist_pt_w_hadronic_zzj_production,
            "hist_delta_phi_lep_met": hist_delta_phi_lep_met_zzj_production,
            "hist_mt_w_leptonic": hist_mt_w_leptonic_zzj_production,
            "hist_ht_total": hist_ht_total_zzj_production,
            "hist_delta_phi_jj": hist_delta_phi_jj_zzj_production,
            "hist_delta_phi_wl_wh": hist_delta_phi_wl_wh_zzj_production,
            "hist_delta_eta_wl_wh": hist_delta_eta_wl_wh_zzj_production,
            "hist_m_jj": hist_m_jj_zzj_production,
            "hist_m_lvjj": hist_m_lvjj_zzj_production
        },
        "wzj": {
            "hist_lepton_pt": hist_lepton_pt_wzj_production,
            "hist_leading_jet_pt": hist_leading_jet_pt_wzj_production,
            "hist_lepton_eta": hist_lepton_eta_wzj_production,
            "hist_delta_r": hist_delta_r_wzj_production,
            "hist_missing_et": hist_missing_et_wzj_production,
            "hist_subleading_jet_eta": hist_subleading_jet_eta_wzj_production,
            "hist_leading_jet_eta": hist_leading_jet_eta_wzj_production,
            "hist_jet_centrality": hist_jet_centrality_wzj_production,
            "hist_delta_eta_jj": hist_delta_eta_jj_wzj_production,
            "hist_m_w_hadronic": hist_m_w_hadronic_wzj_production,
            "hist_m_w_leptonic": hist_m_w_leptonic_wzj_production,
            # ✅ New observables
            "hist_pt_w_leptonic": hist_pt_w_leptonic_wzj_production,
            "hist_pt_w_hadronic": hist_pt_w_hadronic_wzj_production,
            "hist_delta_phi_lep_met": hist_delta_phi_lep_met_wzj_production,
            "hist_mt_w_leptonic": hist_mt_w_leptonic_wzj_production,
            "hist_ht_total": hist_ht_total_wzj_production,
            "hist_delta_phi_jj": hist_delta_phi_jj_wzj_production,
            "hist_delta_phi_wl_wh": hist_delta_phi_wl_wh_wzj_production,
            "hist_delta_eta_wl_wh": hist_delta_eta_wl_wh_wzj_production,
            "hist_m_jj": hist_m_jj_wzj_production,
            "hist_m_lvjj": hist_m_lvjj_wzj_production
        }
    }

    # ✅ Save each background set in its own directory
    for bg_name, histograms in background_histogram_sets.items():
        bg_dir = output_file.mkdir(f"background_{bg_name}")
        bg_dir.cd()
        for hist_name, hist in histograms.items():
            if hist:
                hist.Write()
            else:
                print(f"⚠️ Warning: {hist_name} in {bg_name} is empty!")

    # ✅ Save each signal histogram set
    for signal_name, histograms in signal_histograms.items():
        if signal_name in signal_dirs:
            signal_dirs[signal_name].cd()
            for hist_name, hist in histograms.items():
                if hist:
                    hist.Write()
                else:
                    print(f"⚠️ Warning: {hist_name} for {signal_name} is empty!")

    print("✅ All signal and background histograms written to ROOT file.")

except Exception as e:
    print(f"❌ Error while saving histograms: {e}")

finally:
    output_file.Close()
    print("📁 ROOT file closed successfully.")


#=========================================================================
#=========================================================================



