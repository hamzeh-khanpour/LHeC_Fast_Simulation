

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
    "$FM_{0} / \Lambda^4$": "aa_ww_semi_leptonic_NP_FM0_100K.root",
    "$FM_{1} / \Lambda^4$": "aa_ww_semi_leptonic_NP_FM1_100K.root",
    "$FM_{2} / \Lambda^4$": "aa_ww_semi_leptonic_NP_FM2_100K.root",
    "$FM_{3} / \Lambda^4$": "aa_ww_semi_leptonic_NP_FM3_100K.root",
}
background_file_path = "aa_ww_semi_leptonic_SM_100K.root"




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
        hist_lepton.Fill(leptons[0].Pt())
        hist_lepton_eta.Fill(leptons[0].Eta())   # ✅ Fill η histogram

        # Fill histogram with leading jet pT
        leading_jet = jets[0] if jets[0].Pt() > jets[1].Pt() else jets[1]
        hist_leading_jet.Fill(leading_jet.Pt())


        # Fill optional histograms  Fill histogram with leading jet pT
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
        if leading_jet.Pt() < 100 or leptons[0].Pt() < 100:
            continue

        # Count selected events
        selected_events_final += 1



    # **✅ Selection Efficiency Calculations**
    efficiency_pre = selected_events_pre / total_events if total_events > 0 else 0

    efficiency_final = selected_events_final / total_events if total_events > 0 else 0



    # Return histograms and efficiency
    return (
        hist_lepton,
        hist_leading_jet,
        hist_lepton_eta,
        hist_delta_r,
        hist_missing_et,
        hist_centrality,
        hist_exp_centrality,
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

# Parameters for differential cross-section

signal_cross_sections = {
    "$FM_{0} / \Lambda^4$": 1.99216,   # pb
    "$FM_{1} / \Lambda^4$": 0.17783,   # pb
    "$FM_{2} / \Lambda^4$": 77.8809,   # pb
    "$FM_{3} / \Lambda^4$": 5.95386    # pb
}

background_cross_section = 0.0134936  # pb



#signal_cross_section_0   = 58.8643      # pb  FT0
#signal_cross_section_0   = 6.53478      # pb  FT1
#signal_cross_section_0   = 4.57082      # pb  FT2


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

# Define ranges for invariant masses
m_w_hadronic_range = (1, 150)  # Range for the hadronic W boson mass
m_w_leptonic_range = (1, 150)  # Range for the leptonic W boson mass



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
# Dictionary to store histograms for each signal
signal_histograms = {}

for signal_name in signal_files:
    signal_histograms[signal_name] = {
        "hist_lepton": ROOT.TH1F(f"hist_lepton_{signal_name}", f"Lepton pT Distribution ({signal_name}); p_{{T}} [GeV]; Entries", num_bins, *pt_range_lepton),
        "hist_leading_jet": ROOT.TH1F(f"hist_leading_jet_{signal_name}", f"Leading Jet pT Distribution ({signal_name}); p_{{T}} [GeV]; Entries", num_bins, *pt_range_jet),
        "hist_lepton_eta": ROOT.TH1F(f"hist_lepton_eta_{signal_name}", f"Lepton Eta Distribution ({signal_name}); #eta; Entries", num_bins, *eta_range),
        "hist_delta_r": ROOT.TH1F(f"hist_delta_r_{signal_name}", f"Delta R Distribution ({signal_name}); ΔR; Entries", num_bins, *delta_r_range),
        "hist_missing_et": ROOT.TH1F(f"hist_missing_et_{signal_name}", f"Missing ET Distribution ({signal_name}); MET [GeV]; Entries", num_bins, *met_range),
        "hist_centrality": ROOT.TH1F(f"hist_centrality_{signal_name}", f"Centrality Distribution ({signal_name}); Centrality; Entries", num_bins, *centrality_range),
        "hist_exp_centrality": ROOT.TH1F(f"hist_exp_centrality_{signal_name}", f"Exponential Centrality Distribution ({signal_name}); exp(Centrality); Entries", num_bins, *exp_centrality_range),
        "hist_jet_centrality": ROOT.TH1F(f"hist_jet_centrality_{signal_name}", f"Jet Centrality Distribution ({signal_name}); Jet Centrality; Entries", num_bins, *jet_centrality_range),
        "hist_delta_eta_jj": ROOT.TH1F(f"hist_delta_eta_jj_{signal_name}", f"Delta Eta Between Jets ({signal_name}); Δηjj; Entries", num_bins, *delta_eta_jj_range),
        "hist_m_w_leptonic": ROOT.TH1F(f"hist_m_w_leptonic_{signal_name}", f"Leptonic W Boson Mass ({signal_name}); m_{{W}}^{{leptonic}} [GeV]; Entries", num_bins, *m_w_leptonic_range),
        "hist_m_w_hadronic": ROOT.TH1F(f"hist_m_w_hadronic_{signal_name}", f"Hadronic W Boson Mass ({signal_name}); m_{{W}}^{{hadronic}} [GeV]; Entries", num_bins, *m_w_hadronic_range)
    }

# Background histograms remain the same
hist_lepton_background = ROOT.TH1F("hist_lepton_background", "Lepton pT Distribution (Background); p_{T} [GeV]; Entries", num_bins, *pt_range_lepton)
hist_leading_jet_background = ROOT.TH1F("hist_leading_jet_background", "Leading Jet pT Distribution (Background); p_{T} [GeV]; Entries", num_bins, *pt_range_jet)
hist_lepton_eta_background = ROOT.TH1F("hist_lepton_eta_background", "Lepton Eta Distribution (Background); #eta; Entries", num_bins, *eta_range)
hist_delta_r_background = ROOT.TH1F("hist_delta_r_background", "Delta R Distribution (Background); ΔR; Entries", num_bins, *delta_r_range)
hist_missing_et_background = ROOT.TH1F("hist_missing_et_background", "Missing ET Distribution (Background); MET [GeV]; Entries", num_bins, *met_range)
hist_centrality_background = ROOT.TH1F("hist_centrality_background", "Centrality Distribution (Background); Centrality; Entries", num_bins, *centrality_range)
hist_exp_centrality_background = ROOT.TH1F("hist_exp_centrality_background", "Exponential Centrality Distribution (Background); exp(Centrality); Entries", num_bins, *exp_centrality_range)
hist_jet_centrality_background = ROOT.TH1F("hist_jet_centrality_background", "Jet Centrality Distribution (Background); Jet Centrality; Entries", num_bins, *jet_centrality_range)
hist_delta_eta_jj_background = ROOT.TH1F("hist_delta_eta_jj_background", "Delta Eta Between Jets (Background); Δηjj; Entries", num_bins, *delta_eta_jj_range)
hist_m_w_hadronic_background = ROOT.TH1F("hist_m_w_hadronic_background", "Hadronic W Boson Mass (Background); m_{W}^{hadronic} [GeV]; Entries", num_bins, *m_w_hadronic_range)
hist_m_w_leptonic_background = ROOT.TH1F("hist_m_w_leptonic_background", "Leptonic W Boson Mass (Background); m_{W}^{leptonic} [GeV]; Entries", num_bins, *m_w_leptonic_range)




# Dictionary to store efficiencies for each signal
signal_efficiencies = {}

# Process all signal files dynamically
for signal_name, file_path in signal_files.items():
    print(f"Processing signal: {signal_name}")

    # Get the corresponding histograms for this signal
    histograms = signal_histograms[signal_name]

    # Process the signal file
    (histograms["hist_lepton"], histograms["hist_leading_jet"], histograms["hist_lepton_eta"], histograms["hist_delta_r"],
     histograms["hist_missing_et"], histograms["hist_centrality"], histograms["hist_exp_centrality"], histograms["hist_jet_centrality"],
     histograms["hist_delta_eta_jj"], histograms["hist_m_w_leptonic"], histograms["hist_m_w_hadronic"],
     efficiency_pre, efficiency_final) = process_file(
        file_path, histograms["hist_lepton"], histograms["hist_leading_jet"], histograms["hist_lepton_eta"],
        histograms["hist_delta_r"], histograms["hist_missing_et"], histograms["hist_centrality"], histograms["hist_exp_centrality"],
        histograms["hist_jet_centrality"], histograms["hist_delta_eta_jj"], histograms["hist_m_w_leptonic"], histograms["hist_m_w_hadronic"]
    )

    # Store efficiencies for this signal
    signal_efficiencies[signal_name] = {
        "efficiency_pre": efficiency_pre,
        "efficiency_final": efficiency_final
    }

# Process the background file separately
(hist_lepton_background, hist_leading_jet_background, hist_lepton_eta_background, hist_delta_r_background,
 hist_missing_et_background, hist_centrality_background, hist_exp_centrality_background, hist_jet_centrality_background,
 hist_delta_eta_jj_background, hist_m_w_leptonic_background, hist_m_w_hadronic_background,
 background_efficiency_pre, background_efficiency_final) = process_file(
    background_file_path, hist_lepton_background, hist_leading_jet_background, hist_lepton_eta_background,
    hist_delta_r_background, hist_missing_et_background, hist_centrality_background, hist_exp_centrality_background,
    hist_jet_centrality_background, hist_delta_eta_jj_background, hist_m_w_leptonic_background, hist_m_w_hadronic_background
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




# Dictionary to store differential cross-sections for each signal
signal_dsigma = {}

# Calculate differential cross-sections for all signals
for signal_name, histograms in signal_histograms.items():
    print(f"Calculating differential cross-sections for {signal_name}...")

    cross_section = signal_cross_sections[signal_name]  # Get cross-section for this signal

    signal_dsigma[signal_name] = {
        "pt_bins_lepton": calculate_dsigma(histograms["hist_lepton"], cross_section, bin_width_pt_lepton),
        "pt_bins_jet": calculate_dsigma(histograms["hist_leading_jet"], cross_section, bin_width_pt_jet),
        "eta_bins_lepton": calculate_dsigma(histograms["hist_lepton_eta"], cross_section, bin_width_eta),
        "delta_r_bins": calculate_dsigma(histograms["hist_delta_r"], cross_section, bin_width_delta_r),
        "met_bins": calculate_dsigma(histograms["hist_missing_et"], cross_section, bin_width_met),
        "centrality_bins": calculate_dsigma(histograms["hist_centrality"], cross_section, bin_width_centrality),
        "exp_centrality_bins": calculate_dsigma(histograms["hist_exp_centrality"], cross_section, bin_width_exp_centrality),
        "jet_centrality_bins": calculate_dsigma(histograms["hist_jet_centrality"], cross_section, bin_width_jet_centrality),
        "delta_eta_jj_bins": calculate_dsigma(histograms["hist_delta_eta_jj"], cross_section, bin_width_delta_eta_jj),
        "m_w_hadronic_bins": calculate_dsigma(histograms["hist_m_w_hadronic"], cross_section, bin_width_m_w_hadronic),
        "m_w_leptonic_bins": calculate_dsigma(histograms["hist_m_w_leptonic"], cross_section, bin_width_m_w_leptonic),
    }



# Calculate differential cross-sections for background
background_dsigma = {
    "pt_bins_lepton": calculate_dsigma(hist_lepton_background, background_cross_section, bin_width_pt_lepton),
    "pt_bins_jet": calculate_dsigma(hist_leading_jet_background, background_cross_section, bin_width_pt_jet),
    "eta_bins_lepton": calculate_dsigma(hist_lepton_eta_background, background_cross_section, bin_width_eta),
    "delta_r_bins": calculate_dsigma(hist_delta_r_background, background_cross_section, bin_width_delta_r),
    "met_bins": calculate_dsigma(hist_missing_et_background, background_cross_section, bin_width_met),
    "centrality_bins": calculate_dsigma(hist_centrality_background, background_cross_section, bin_width_centrality),
    "exp_centrality_bins": calculate_dsigma(hist_exp_centrality_background, background_cross_section, bin_width_exp_centrality),
    "jet_centrality_bins": calculate_dsigma(hist_jet_centrality_background, background_cross_section, bin_width_jet_centrality),
    "delta_eta_jj_bins": calculate_dsigma(hist_delta_eta_jj_background, background_cross_section, bin_width_delta_eta_jj),
    "m_w_hadronic_bins": calculate_dsigma(hist_m_w_hadronic_background, background_cross_section, bin_width_m_w_hadronic),
    "m_w_leptonic_bins": calculate_dsigma(hist_m_w_leptonic_background, background_cross_section, bin_width_m_w_leptonic),
}







#=========================================================================
#=========================================================================




plt.figure(figsize=(11, 12))  # Create a new figure for the leading jet η plot
plt.subplots_adjust(left=0.15, right=0.95, bottom=0.12, top=0.95)



# Define colors for each signal
signal_colors = {
    "$FM_{0} / \Lambda^4$": "red",
    "$FM_{1} / \Lambda^4$": "purple",
    "$FM_{2} / \Lambda^4$": "green",
    "$FM_{3} / \Lambda^4$": "orange"
}


# Create figure for lepton pT differential cross-section

# Loop through all signals and plot their differential cross-section
for signal_name, dsigma_data in signal_dsigma.items():
    pt_bins, dsigma = dsigma_data["pt_bins_lepton"]
    plt.step(pt_bins, dsigma, where="mid", alpha=0.7,
             label=f"LHeC@1.2 TeV : Signal ($w^+ w^-$) [{signal_name}]",
             color=signal_colors.get(signal_name, "black"), linewidth=3)

# Plot the background as well
pt_bins_background, dsigma_background = background_dsigma["pt_bins_lepton"]
plt.step(pt_bins_background, dsigma_background, where="mid", alpha=0.7,
         label="LHeC@1.2 TeV : SM background ($w^+ w^-$)", color="blue", linewidth=3)

# Set labels and title
plt.xlabel(r"$p_T^{\ell} \ \mathrm{[GeV]}$")
plt.ylabel(r"$\frac{d\sigma}{dp_T^{\ell}} \ \mathrm{[pb/GeV]}$")
plt.title(r"Delphes simulation : $e^- p \to e^- w^+ w^- p \to e^- j j \ell \nu_{\ell} p$ : LHeC@1.2 TeV", fontsize=20)
plt.yscale("log")

# Grid, legend, and layout adjustments
plt.legend()
plt.grid(True, linestyle="--", alpha=0.6)
plt.tight_layout()
plt.ylim(0.00001, 10.0)

# Save the plot
plt.savefig("/home/hamzeh-khanpour/Documents/GitHub/LHeC_Fast_Simulation/Pythia8_Delphes/differential_cross_section_lepton_pt_allsignal.pdf", dpi=600)

# Show the plot
plt.show()







# Create figure for leading jet pT differential cross-section
plt.figure(figsize=(11, 12))  # Create a new figure for the leading jet η plot

# Loop through all signals and plot their differential cross-section
for signal_name, dsigma_data in signal_dsigma.items():
    pt_bins, dsigma = dsigma_data["pt_bins_jet"]
    plt.step(pt_bins, dsigma, where="mid", alpha=0.7,
             label=f"LHeC@1.2 TeV : Signal ($w^+ w^-$) [{signal_name}]",
             color=signal_colors.get(signal_name, "black"), linewidth=3)

# Plot the background as well
pt_bins_background, dsigma_background = background_dsigma["pt_bins_jet"]
plt.step(pt_bins_background, dsigma_background, where="mid", alpha=0.7,
         label="LHeC@1.2 TeV : SM background ($w^+ w^-$)", color="blue", linewidth=3)

# Set labels and title
plt.xlabel(r"$p_T^{\mathrm{leading~jet}} \ \mathrm{[GeV]}$")
plt.ylabel(r"$\frac{d\sigma}{dp_T^{\mathrm{leading~jet}}} \ \mathrm{[pb/GeV]}$")
plt.title(r"Delphes simulation : $e^- p \to e^- w^+ w^- p \to e^- j j \ell \nu_{\ell} p$ : LHeC@1.2 TeV", fontsize=20)
plt.yscale("log")

# Grid, legend, and layout adjustments
plt.legend()
plt.grid(True, linestyle="--", alpha=0.6)
plt.tight_layout()
plt.ylim(0.00001, 10.0)

# Save the plot
plt.savefig("/home/hamzeh-khanpour/Documents/GitHub/LHeC_Fast_Simulation/Pythia8_Delphes/differential_cross_section_jet_pt_allsignal.pdf", dpi=600)

# Show the plot
plt.show()









# Plot the differential cross-sections for lepton η
plt.figure(figsize=(11, 12))  # Create a new figure for the leading jet η plot


# Loop through all signals and plot their differential cross-section
for signal_name, dsigma_data in signal_dsigma.items():
    eta_bins, dsigma = dsigma_data["eta_bins_lepton"]
    plt.step(eta_bins, dsigma, where="mid", alpha=0.7,
             label=f"LHeC@1.2 TeV : Signal ($w^+ w^-$) [{signal_name}]",
             color=signal_colors.get(signal_name, "black"), linewidth=3)

# Plot the background as well
eta_bins_background, dsigma_background = background_dsigma["eta_bins_lepton"]
plt.step(eta_bins_background, dsigma_background, where="mid", alpha=0.7,
         label="LHeC@1.2 TeV : SM background ($w^+ w^-$)", color="blue", linewidth=3)

# Set labels and title
plt.xlabel(r"$\eta^{\ell}$")
plt.ylabel(r"$\frac{d\sigma}{d\eta^{\ell}} \ \mathrm{[pb]}$")
plt.title(r"Delphes simulation : $e^- p \to e^- w^+ w^- p \to e^- j j \ell \nu_{\ell} p$ : LHeC@1.2 TeV", fontsize=20)
plt.yscale("log")

# Grid, legend, and layout adjustments
plt.legend()
plt.grid(True, linestyle="--", alpha=0.6)
plt.tight_layout()
plt.ylim(0.0001, 1000.0)

# Save the plot
plt.savefig("/home/hamzeh-khanpour/Documents/GitHub/LHeC_Fast_Simulation/Pythia8_Delphes/differential_cross_section_lepton_eta_allsignal.pdf", dpi=600)

# Show the plot
plt.show()








# Create figure for Delta R differential cross-section
plt.figure(figsize=(11, 12))  # Create a new figure for the leading jet η plot

# Loop through all signals and plot their differential cross-section
for signal_name, dsigma_data in signal_dsigma.items():
    delta_r_bins, dsigma = dsigma_data["delta_r_bins"]
    plt.step(delta_r_bins, dsigma, where="mid", alpha=0.7,
             label=f"LHeC@1.2 TeV : Signal ($w^+ w^-$) [{signal_name}]",
             color=signal_colors.get(signal_name, "black"), linewidth=3)

# Plot the background as well
delta_r_bins_background, dsigma_background = background_dsigma["delta_r_bins"]
plt.step(delta_r_bins_background, dsigma_background, where="mid", alpha=0.7,
         label="LHeC@1.2 TeV : SM background ($w^+ w^-$)", color="blue", linewidth=3)

# Set labels and title
plt.xlabel(r"$\Delta R(\ell, \mathrm{leading~jet})$")
plt.ylabel(r"$\frac{d\sigma}{d\Delta R} \ \mathrm{[pb]}$")
plt.title(r"Delphes simulation : $e^- p \to e^- w^+ w^- p \to e^- j j \ell \nu_{\ell} p$ : LHeC@1.2 TeV", fontsize=20)
plt.yscale("log")

# Grid, legend, and layout adjustments
plt.legend()
plt.grid(True, linestyle="--", alpha=0.6)
plt.tight_layout()
plt.ylim(0.0001, 1000.0)

# Save the plot
plt.savefig("/home/hamzeh-khanpour/Documents/GitHub/LHeC_Fast_Simulation/Pythia8_Delphes/differential_cross_section_delta_r_allsignal.pdf", dpi=600)

# Show the plot
plt.show()








# Create figure for MET differential cross-section
plt.figure(figsize=(11, 12))  # Create a new figure for the leading jet η plot

# Loop through all signals and plot their differential cross-section
for signal_name, dsigma_data in signal_dsigma.items():
    met_bins, dsigma = dsigma_data["met_bins"]
    plt.step(met_bins, dsigma, where="mid", alpha=0.7,
             label=f"LHeC@1.2 TeV : Signal ($w^+ w^-$) [{signal_name}]",
             color=signal_colors.get(signal_name, "black"), linewidth=3)

# Plot the background as well
met_bins_background, dsigma_background = background_dsigma["met_bins"]
plt.step(met_bins_background, dsigma_background, where="mid", alpha=0.7,
         label="LHeC@1.2 TeV : SM background ($w^+ w^-$)", color="blue", linewidth=3)

# Set labels and title
plt.xlabel(r"$\mathrm{MET} \ \mathrm{[GeV]}$")
plt.ylabel(r"$\frac{d\sigma}{d\mathrm{MET}} \ \mathrm{[pb/GeV]}$")
plt.title(r"Delphes simulation : $e^- p \to e^- w^+ w^- p \to e^- j j \ell \nu_{\ell} p$ : LHeC@1.2 TeV", fontsize=20)
plt.yscale("log")

# Grid, legend, and layout adjustments
plt.legend()
plt.grid(True, linestyle="--", alpha=0.6)
plt.tight_layout()
plt.ylim(0.00001, 10.0)

# Save the plot
plt.savefig("/home/hamzeh-khanpour/Documents/GitHub/LHeC_Fast_Simulation/Pythia8_Delphes/differential_cross_section_met_allsignal.pdf", dpi=600)

# Show the plot
plt.show()







# Create figure for normalized hadronic W boson mass
plt.figure(figsize=(11, 12))  # Create a new figure for the leading jet η plot

# Loop through all signals and normalize their distributions
for signal_name, dsigma_data in signal_dsigma.items():
    m_w_hadronic_bins, dsigma = dsigma_data["m_w_hadronic_bins"]
    normalized_dsigma = dsigma / np.sum(dsigma)  # Normalize histogram
    plt.step(m_w_hadronic_bins, normalized_dsigma, where="mid", alpha=0.7,
             label=f"LHeC@1.2 TeV : Signal ($w^+ w^-$) [{signal_name}]",
             color=signal_colors.get(signal_name, "black"), linewidth=3)

# Plot the background as well (normalized)
m_w_hadronic_bins_background, dsigma_background = background_dsigma["m_w_hadronic_bins"]
normalized_dsigma_background = dsigma_background / np.sum(dsigma_background)
plt.step(m_w_hadronic_bins_background, normalized_dsigma_background, where="mid", alpha=0.7,
         label="LHeC@1.2 TeV : SM background ($w^+ w^-$)", color="blue", linewidth=3)

# Set labels and title
plt.xlabel(r"$M_W^{\mathrm{j_1j_2}} \ \mathrm{[GeV]}$")
plt.ylabel("Normalized Entries")
plt.title(r"Delphes simulation : $e^- p \to e^- w^+ w^- p \to e^- j j \ell \nu_{\ell} p$ : LHeC@1.2 TeV", fontsize=20)

# Grid, legend, and layout adjustments
plt.legend()
plt.grid(True, linestyle="--", alpha=0.6)
plt.tight_layout()
plt.ylim(0.0, 0.3)

# Save the plot
plt.savefig("/home/hamzeh-khanpour/Documents/GitHub/LHeC_Fast_Simulation/Pythia8_Delphes/normalized_m_w_hadronic_allsignal.pdf", dpi=600)

# Show the plot
plt.show()











# Create figure for normalized leptonic W boson mass
plt.figure(figsize=(11, 12))  # Create a new figure for the leading jet η plot

# Loop through all signals and normalize their distributions
for signal_name, dsigma_data in signal_dsigma.items():
    m_w_leptonic_bins, dsigma = dsigma_data["m_w_leptonic_bins"]
    normalized_dsigma = dsigma / np.sum(dsigma)  # Normalize histogram
    plt.step(m_w_leptonic_bins, normalized_dsigma, where="mid", alpha=0.7,
             label=f"LHeC@1.2 TeV : Signal ($w^+ w^-$) [{signal_name}]",
             color=signal_colors.get(signal_name, "black"), linewidth=3)

# Plot the background as well (normalized)
m_w_leptonic_bins_background, dsigma_background = background_dsigma["m_w_leptonic_bins"]
normalized_dsigma_background = dsigma_background / np.sum(dsigma_background)
plt.step(m_w_leptonic_bins_background, normalized_dsigma_background, where="mid", alpha=0.7,
         label="LHeC@1.2 TeV : SM background ($w^+ w^-$)", color="blue", linewidth=3)

# Set labels and title
plt.xlabel(r"$M_W^{\ell \nu_\ell} \ \mathrm{[GeV]}$")
plt.ylabel("Normalized Entries")
plt.title(r"Delphes simulation : $e^- p \to e^- w^+ w^- p \to e^- j j \ell \nu_{\ell} p$ : LHeC@1.2 TeV", fontsize=20)

# Grid, legend, and layout adjustments
plt.legend()
plt.grid(True, linestyle="--", alpha=0.6)
plt.tight_layout()
plt.ylim(0.0, 0.3)

# Save the plot
plt.savefig("/home/hamzeh-khanpour/Documents/GitHub/LHeC_Fast_Simulation/Pythia8_Delphes/normalized_m_w_leptonic_allsignal.pdf", dpi=600)

# Show the plot
plt.show()







#=========================================================================
#=========================================================================

# Open the output ROOT file (use "UPDATE" if you want to append)
output_file = ROOT.TFile("/home/hamzeh-khanpour/Documents/GitHub/LHeC_Fast_Simulation/Pythia8_Delphes/output_histograms.root", "RECREATE")

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
        "hist_centrality": hist_centrality_background,
        "hist_exp_centrality": hist_exp_centrality_background,
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



