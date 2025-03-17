
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
signal_file_path = "aa_ww_semi_leptonic_NP_1_FM2_100.root"
background_file_path = "aa_ww_semi_leptonic_SM.root"



# Load Delphes library
ROOT.gSystem.Load("libDelphes")
ROOT.gInterpreter.Declare('#include "classes/DelphesClasses.h"')
ROOT.gInterpreter.Declare('#include "external/ExRootAnalysis/ExRootTreeReader.h"')


## ================================================================================



# Function to process a ROOT file and fill histograms for leptons and leading jet
def process_file(
    file_path,
    hist_lepton,
    hist_leading_jet,
    hist_m_w_leptonic,
    hist_m_w_hadronic,
    hist_lepton_eta,
    hist_leading_jet_eta,
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

        # Fill histogram with lepton pT
        hist_lepton.Fill(leptons[0].Pt())
        hist_lepton_eta.Fill(leptons[0].Eta())  # ✅ Fill η histogram

        # Fill histogram with leading jet pT
        leading_jet = max(jets, key=lambda jet: jet.Pt())  # Select highest pT jet
        hist_leading_jet.Fill(leading_jet.Pt())

        hist_leading_jet_eta.Fill(leading_jet.Eta())



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


        # **✅ Hadronic W Reconstruction**
        if hist_m_w_hadronic is not None:
            w_hadronic = jets[0] + jets[1]
            hist_m_w_hadronic.Fill(w_hadronic.M())


    # Calculate selection efficiency
    efficiency = selected_events / total_events if total_events > 0 else 0

    return hist_lepton, hist_leading_jet, hist_m_w_leptonic, hist_m_w_hadronic, hist_lepton_eta, hist_leading_jet_eta, efficiency




## ================================================================================




# Create histograms for lepton and leading jet pT (signal and background)
hist_lepton_signal = ROOT.TH1F("hist_lepton_signal", "Lepton pT Distribution; p_{T} [GeV]; Entries", 50, 0, 300)
hist_lepton_background = ROOT.TH1F("hist_lepton_background", "Lepton pT Distribution; p_{T} [GeV]; Entries", 50, 0, 300)


hist_leading_jet_signal = ROOT.TH1F("hist_leading_jet_signal", "Leading Jet pT Distribution; p_{T} [GeV]; Entries", 50, 0, 300)
hist_leading_jet_background = ROOT.TH1F("hist_leading_jet_background", "Leading Jet pT Distribution; p_{T} [GeV]; Entries", 50, 0, 300)


hist_m_w_leptonic_signal = ROOT.TH1F("hist_m_w_leptonic_signal", "Leptonic W Mass Distribution (Signal); M_W^{\ell\nu} [GeV]; Entries", 50, 1, 130)
hist_m_w_leptonic_background = ROOT.TH1F("hist_m_w_leptonic_background", "Leptonic W Mass Distribution (Background); M_W^{\ell\nu} [GeV]; Entries", 50, 1, 130)


hist_m_w_hadronic_signal = ROOT.TH1F("hist_m_w_hadronic_signal", "Hadronic W Mass Distribution (Signal); M_W^{jj} [GeV]; Entries", 50, 1, 130)
hist_m_w_hadronic_background = ROOT.TH1F("hist_m_w_hadronic_background", "Hadronic W Mass Distribution (Background); M_W^{jj} [GeV]; Entries", 50, 1, 130)


# ✅ NEW: Define histograms for lepton pseudorapidity (η)
hist_lepton_eta_signal = ROOT.TH1F("hist_lepton_eta_signal", "Lepton Pseudorapidity Distribution (Signal); \eta_{\ell}; Entries", 50, -5, 5)
hist_lepton_eta_background = ROOT.TH1F("hist_lepton_eta_background", "Lepton Pseudorapidity Distribution (Background); \eta_{\ell}; Entries", 50, -5, 5)


# ✅ NEW: Define histograms for leading jet pseudorapidity (η)
hist_leading_jet_eta_signal = ROOT.TH1F("hist_leading_jet_eta_signal", "Leading Jet Pseudorapidity Distribution (Signal); \eta_{\mathrm{jet}}; Entries", 50, -5, 5)
hist_leading_jet_eta_background = ROOT.TH1F("hist_leading_jet_eta_background", "Leading Jet Pseudorapidity Distribution (Background); \eta_{\mathrm{jet}}; Entries", 50, -5, 5)






# Process signal and background files and calculate efficiencies
# Call the process_file function and update the unpacking to include all returned values

(hist_lepton_signal,
 hist_leading_jet_signal,
 hist_m_w_leptonic_signal,
 hist_m_w_hadronic_signal,
 hist_lepton_eta_signal,   # ✅ NEW: Lepton η histogram (Signal)
 hist_leading_jet_eta_signal,  # ✅ NEW: Leading Jet η histogram (Signal)
 signal_efficiency) = process_file(
    signal_file_path,
    hist_lepton_signal,
    hist_leading_jet_signal,
    hist_m_w_leptonic_signal,
    hist_m_w_hadronic_signal,
    hist_lepton_eta_signal,   # ✅ Pass to function
    hist_leading_jet_eta_signal  # ✅ Pass to function
)


(hist_lepton_background,
 hist_leading_jet_background,
 hist_m_w_leptonic_background,
 hist_m_w_hadronic_background,
 hist_lepton_eta_background,   # ✅ NEW: Lepton η histogram (Background)
 hist_leading_jet_eta_background,  # ✅ NEW: Leading Jet η histogram (Background)
 background_efficiency) = process_file(
    background_file_path,
    hist_lepton_background,
    hist_leading_jet_background,
    hist_m_w_leptonic_background,
    hist_m_w_hadronic_background,
    hist_lepton_eta_background,   # ✅ Pass to function
    hist_leading_jet_eta_background  # ✅ Pass to function
)





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


# Convert ROOT histograms to matplotlib for leptonic W boson invariant mass
x_vals_m_w_leptonic_signal = [hist_m_w_leptonic_signal.GetBinCenter(i) for i in range(1, hist_m_w_leptonic_signal.GetNbinsX() + 1)]
y_vals_m_w_leptonic_signal = [hist_m_w_leptonic_signal.GetBinContent(i) for i in range(1, hist_m_w_leptonic_signal.GetNbinsX() + 1)]

x_vals_m_w_leptonic_background = [hist_m_w_leptonic_background.GetBinCenter(i) for i in range(1, hist_m_w_leptonic_background.GetNbinsX() + 1)]
y_vals_m_w_leptonic_background = [hist_m_w_leptonic_background.GetBinContent(i) for i in range(1, hist_m_w_leptonic_background.GetNbinsX() + 1)]



# Convert ROOT histograms to matplotlib for hadronic W boson invariant mass
x_vals_m_w_hadronic_signal = [hist_m_w_hadronic_signal.GetBinCenter(i) for i in range(1, hist_m_w_hadronic_signal.GetNbinsX() + 1)]
y_vals_m_w_hadronic_signal = [hist_m_w_hadronic_signal.GetBinContent(i) for i in range(1, hist_m_w_hadronic_signal.GetNbinsX() + 1)]

x_vals_m_w_hadronic_background = [hist_m_w_hadronic_background.GetBinCenter(i) for i in range(1, hist_m_w_hadronic_background.GetNbinsX() + 1)]
y_vals_m_w_hadronic_background = [hist_m_w_hadronic_background.GetBinContent(i) for i in range(1, hist_m_w_hadronic_background.GetNbinsX() + 1)]




# ✅ NEW: Convert ROOT histograms for lepton pseudorapidity (η)
x_vals_lepton_eta_signal = [hist_lepton_eta_signal.GetBinCenter(i) for i in range(1, hist_lepton_eta_signal.GetNbinsX() + 1)]
y_vals_lepton_eta_signal = [hist_lepton_eta_signal.GetBinContent(i) for i in range(1, hist_lepton_eta_signal.GetNbinsX() + 1)]

x_vals_lepton_eta_background = [hist_lepton_eta_background.GetBinCenter(i) for i in range(1, hist_lepton_eta_background.GetNbinsX() + 1)]
y_vals_lepton_eta_background = [hist_lepton_eta_background.GetBinContent(i) for i in range(1, hist_lepton_eta_background.GetNbinsX() + 1)]



# ✅ NEW: Convert ROOT histograms for leading jet pseudorapidity (η)
x_vals_leading_jet_eta_signal = [hist_leading_jet_eta_signal.GetBinCenter(i) for i in range(1, hist_leading_jet_eta_signal.GetNbinsX() + 1)]
y_vals_leading_jet_eta_signal = [hist_leading_jet_eta_signal.GetBinContent(i) for i in range(1, hist_leading_jet_eta_signal.GetNbinsX() + 1)]

x_vals_leading_jet_eta_background = [hist_leading_jet_eta_background.GetBinCenter(i) for i in range(1, hist_leading_jet_eta_background.GetNbinsX() + 1)]
y_vals_leading_jet_eta_background = [hist_leading_jet_eta_background.GetBinContent(i) for i in range(1, hist_leading_jet_eta_background.GetNbinsX() + 1)]





# Plot the histograms for leptons
plt.figure(figsize=(10, 12))  # Create a new figure for the lepton plot
plt.subplots_adjust(left=0.15, right=0.95, bottom=0.12, top=0.95)

plt.bar(    x_vals_lepton_signal,    y_vals_lepton_signal,    width=hist_lepton_signal.GetBinWidth(1),    alpha=0.6,
    label="LHeC@1.2 TeV : Signal ($w^+ w^-) [f_{M_2} / \Lambda^4$)", color="red",)
plt.bar(    x_vals_lepton_background,    y_vals_lepton_background,    width=hist_lepton_background.GetBinWidth(1),    alpha=0.6,
    label="LHeC@1.2 TeV : SM background ($w^+ w^-$)", color="blue",)
plt.xlabel(r"$p_T^{\ell} \ \mathrm{[GeV]}$")
plt.ylabel("Entries")
plt.title(r"Delphes simulation : $e^- p \to e^- w^+ w^- p \to e^- j j \ell \nu_{\ell} p$ : LHeC@1.2 TeV", fontsize=18)
plt.legend()
plt.grid()
plt.savefig("/home/hamzeh-khanpour/Documents/GitHub/LHeC_Fast_Simulation/Pythia8_Delphes_semi_leptonic/lepton_pt_comparison.pdf", dpi=300)
plt.show()  # Ensure this ends the current plot properly





# Plot the histograms for leading jets
plt.figure(figsize=(10, 12))  # Create a new figure for the jet plot
plt.subplots_adjust(left=0.15, right=0.95, bottom=0.12, top=0.95)

plt.bar(    x_vals_jet_signal,    y_vals_jet_signal,    width=hist_leading_jet_signal.GetBinWidth(1),    alpha=0.6,    label="LHeC@1.2 TeV : Signal ($w^+ w^-) [f_{M_2} / \Lambda^4$)", color="red",)
plt.bar(    x_vals_jet_background,    y_vals_jet_background,    width=hist_leading_jet_background.GetBinWidth(1),    alpha=0.6,
    label="LHeC@1.2 TeV : SM background ($w^+ w^-$)", color="blue",)
plt.xlabel(r"$p_T^{\mathrm{leading~jet}} \ \mathrm{[GeV]}$")
plt.ylabel("Entries")
plt.title(r"Delphes simulation : $e^- p \to e^- w^+ w^- p \to e^- j j \ell \nu_{\ell} p$ : LHeC@1.2 TeV", fontsize=18)
plt.legend()
plt.grid()
plt.savefig("/home/hamzeh-khanpour/Documents/GitHub/LHeC_Fast_Simulation/Pythia8_Delphes_semi_leptonic/leading_jet_pt_comparison.pdf", dpi=300)
plt.show()  # Ensure the second plot is displayed






# Plot the histograms for leptonic W boson invariant mass
plt.figure(figsize=(10, 12))  # Create a new figure for leptonic W mass plot
plt.subplots_adjust(left=0.15, right=0.95, bottom=0.12, top=0.95)

plt.bar(    x_vals_m_w_leptonic_signal,    y_vals_m_w_leptonic_signal,    width=hist_m_w_leptonic_signal.GetBinWidth(1), alpha=0.6, label="LHeC@1.2 TeV : Signal ($w^+ w^-) [f_{M_2} / \Lambda^4$)", color="red",)
plt.bar(    x_vals_m_w_leptonic_background,    y_vals_m_w_leptonic_background,    width=hist_m_w_leptonic_background.GetBinWidth(1), alpha=0.6, label="LHeC@1.2 TeV : SM background ($w^+ w^-$)", color="blue",)
plt.xlabel(r"$M_W^{\ell\nu_{\ell}} \ \mathrm{[GeV]}$")
plt.ylabel("Entries")
plt.title(r"Delphes simulation : $e^- p \to e^- w^+ w^- p \to e^- j j \ell \nu_{\ell} p$ : LHeC@1.2 TeV", fontsize=18)
plt.legend()
plt.grid()
plt.savefig("/home/hamzeh-khanpour/Documents/GitHub/LHeC_Fast_Simulation/Pythia8_Delphes_semi_leptonic/m_w_leptonic_comparison.pdf", dpi=300)
plt.show()





# Plot the histograms for hadronic W boson invariant mass
plt.figure(figsize=(10, 12))  # Create a new figure for hadronic W mass plot
plt.subplots_adjust(left=0.15, right=0.95, bottom=0.12, top=0.95)

plt.bar(    x_vals_m_w_hadronic_signal,    y_vals_m_w_hadronic_signal,    width=hist_m_w_hadronic_signal.GetBinWidth(1), alpha=0.6, label="LHeC@1.2 TeV : Signal ($w^+ w^-) [f_{M_2} / \Lambda^4$)", color="red",)
plt.bar(    x_vals_m_w_hadronic_background,    y_vals_m_w_hadronic_background,    width=hist_m_w_hadronic_background.GetBinWidth(1), alpha=0.6, label="LHeC@1.2 TeV : SM background ($w^+ w^-$)", color="blue",)
plt.xlabel(r"$M_W^{\mathrm{j_1j_2}} \ \mathrm{[GeV]}$")
plt.ylabel("Entries")
plt.title(r"Delphes simulation : $e^- p \to e^- w^+ w^- p \to e^- j j \ell \nu_{\ell} p$ : LHeC@1.2 TeV", fontsize=18)
plt.legend()
plt.grid()
plt.savefig("/home/hamzeh-khanpour/Documents/GitHub/LHeC_Fast_Simulation/Pythia8_Delphes_semi_leptonic/m_w_hadronic_comparison.pdf", dpi=300)
plt.show()





# ✅ NEW: Plot the histograms for lepton pseudorapidity (η)
plt.figure(figsize=(10, 12))  # Create a new figure for the lepton η plot
plt.subplots_adjust(left=0.15, right=0.95, bottom=0.12, top=0.95)

plt.bar(x_vals_lepton_eta_signal, y_vals_lepton_eta_signal, width=hist_lepton_eta_signal.GetBinWidth(1), alpha=0.6, label="LHeC@1.2 TeV : Signal ($w^+ w^-) [f_{M_2} / \Lambda^4$)", color="red")
plt.bar(x_vals_lepton_eta_background, y_vals_lepton_eta_background, width=hist_lepton_eta_background.GetBinWidth(1), alpha=0.6, label="LHeC@1.2 TeV : SM background ($w^+ w^-$)", color="blue")
plt.xlabel(r"$\eta_{\ell}$")
plt.ylabel("Entries")
plt.title(r"Lepton Pseudorapidity Distribution: $e^- p \to e^- w^+ w^- p$", fontsize=18)
plt.legend()
plt.grid()
plt.savefig("/home/hamzeh-khanpour/Documents/GitHub/LHeC_Fast_Simulation/Pythia8_Delphes_semi_leptonic/lepton_eta_comparison.pdf", dpi=300)
plt.show()






# ✅ NEW: Plot the histograms for leading jet pseudorapidity (η)
plt.figure(figsize=(10, 12))  # Create a new figure for the leading jet η plot
plt.subplots_adjust(left=0.15, right=0.95, bottom=0.12, top=0.95)

plt.bar(x_vals_leading_jet_eta_signal, y_vals_leading_jet_eta_signal, width=hist_leading_jet_eta_signal.GetBinWidth(1), alpha=0.6, label="LHeC@1.2 TeV : Signal ($w^+ w^-) [f_{M_2} / \Lambda^4$)", color="red")
plt.bar(x_vals_leading_jet_eta_background, y_vals_leading_jet_eta_background, width=hist_leading_jet_eta_background.GetBinWidth(1), alpha=0.6, label="LHeC@1.2 TeV : SM background ($w^+ w^-$)", color="blue")
plt.xlabel(r"$\eta_{\mathrm{leading~jet}}$")
plt.ylabel("Entries")
plt.title(r"Leading Jet Pseudorapidity Distribution: $e^- p \to e^- w^+ w^- p$", fontsize=18)
plt.legend()
plt.grid()
plt.savefig("/home/hamzeh-khanpour/Documents/GitHub/LHeC_Fast_Simulation/Pythia8_Delphes_semi_leptonic/leading_jet_eta_comparison.pdf", dpi=300)
plt.show()








