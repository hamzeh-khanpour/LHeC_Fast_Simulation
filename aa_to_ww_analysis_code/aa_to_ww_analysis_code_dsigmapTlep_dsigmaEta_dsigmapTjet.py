
# Hamzeh khanpour -- 2025

import ROOT
import matplotlib.pyplot as plt
import numpy as np


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



# Function to parse the LHE file and extract lepton transverse momentum (pT) and leading jet pT
def parse_lhe_file(file_name):
    pt_leptons = []
    eta_leptons = []
    pt_leading_jet = []  # To store the transverse momentum of the leading jet

    with open(file_name, "r") as file:
        in_event = False

        for line in file:
            line = line.strip()

            # Check if we are inside an <event> block
            if "<event>" in line:
                in_event = True
                jets = []  # Temporary list to store jet transverse momenta in an event
                continue
            if "</event>" in line:
                in_event = False
                # Determine the leading jet in this event (if any)
                if jets:
                    pt_leading_jet.append(max(jets))
                continue

            if in_event:
                # Skip lines starting with non-numeric characters (metadata lines)
                if not line[0].isdigit() and not line[0] == "-":  # Include negative PDG IDs
                    continue

                # Split the line into components
                parts = line.split()
                # Check if the line has enough elements to extract particle data
                if len(parts) < 10:
                    continue

                pdg_id = int(parts[0])
                px = float(parts[6])
                py = float(parts[7])
                pz = float(parts[8])
                energy = float(parts[9])

                # Check if the particle is a lepton (e- or e+)
                if abs(pdg_id) in [11, 13, 15]:  # e, mu, tau
                    # Create TLorentzVector
                    lepton = ROOT.TLorentzVector()
                    lepton.SetPxPyPzE(px, py, pz, energy)

                    # Extract transverse momentum
                    pt_leptons.append(lepton.Pt())
                    eta_leptons.append(lepton.Eta())

                # Check if the particle is a jet (status 1, generic criterion for outgoing particles)
                if abs(pdg_id) not in [11, 12, 13, 14, 15, 16, 22] and int(parts[1]) == 1:
                    # Create TLorentzVector for jet
                    jet = ROOT.TLorentzVector()
                    jet.SetPxPyPzE(px, py, pz, energy)

                    # Extract transverse momentum
                    jets.append(jet.Pt())

    return pt_leptons, eta_leptons, pt_leading_jet


#fig, ax = plt.subplots(figsize=(12.0, 10.0))
plt.subplots_adjust(left=0.15, right=0.95, bottom=0.12, top=0.95)



# Parse signal and background files
signal_file_0 = "/home/hamzeh-khanpour/MG5_aMC_v3_5_7/aa_ww_semi_leptonic_NP_FM0/Events/run_01/aa_ww_semi_leptonic_NP_FM0.lhe"
signal_file_2 = "/home/hamzeh-khanpour/MG5_aMC_v3_5_7/aa_ww_semi_leptonic_NP_FM2/Events/run_01/aa_ww_semi_leptonic_NP_FM2.lhe"
background_file = "/home/hamzeh-khanpour/MG5_aMC_v3_5_7/aa_ww_semi_leptonic_SM/Events/run_01/aa_ww_semi_leptonic_SM.lhe"


pt_leptons_signal_0, eta_leptons_signal_0, pt_leading_jet_signal_0 = parse_lhe_file(signal_file_0)
pt_leptons_signal_2, eta_leptons_signal_2, pt_leading_jet_signal_2 = parse_lhe_file(signal_file_2)

pt_leptons_background, eta_leptons_background, pt_leading_jet_background = parse_lhe_file(background_file)



# Parameters for differential cross-section
signal_cross_section_0 = 1.994  # pb
signal_cross_section_2 = 77.86  # pb
background_cross_section = 0.0134802  # pb
num_bins = 50
pt_range_lepton = (0, 400)  # Range for lepton pT
pt_range_jet = (0, 600)     # Range for leading jet pT (adjusted for higher jet momenta)
eta_range = (-10, 10)       # Range for pseudorapidity



# Calculate bin width
bin_width_pt_lepton = (pt_range_lepton[1] - pt_range_lepton[0]) / num_bins
bin_width_pt_jet = (pt_range_jet[1] - pt_range_jet[0]) / num_bins
bin_width_eta = (eta_range[1] - eta_range[0]) / num_bins


# Normalize histograms to calculate differential cross-section
def calculate_dsigma(data, total_cross_section, bin_width, data_range):
    counts, bin_edges = np.histogram(data, bins=num_bins, range=data_range)
    dsigma = counts * (total_cross_section / len(data)) / bin_width
    return bin_edges[:-1], dsigma



# Calculate differential cross-sections for lepton pT, eta, and leading jet pT
pt_bins_signal_0, dsigma_signal_pt_0 = calculate_dsigma(pt_leptons_signal_0, signal_cross_section_0, bin_width_pt_lepton, pt_range_lepton)
pt_bins_signal_2, dsigma_signal_pt_2 = calculate_dsigma(pt_leptons_signal_2, signal_cross_section_2, bin_width_pt_lepton, pt_range_lepton)
pt_bins_background, dsigma_background_pt = calculate_dsigma(pt_leptons_background, background_cross_section, bin_width_pt_lepton, pt_range_lepton)


eta_bins_signal_0, dsigma_signal_eta_0 = calculate_dsigma(eta_leptons_signal_0, signal_cross_section_0, bin_width_eta, eta_range)
eta_bins_signal_2, dsigma_signal_eta_2 = calculate_dsigma(eta_leptons_signal_2, signal_cross_section_2, bin_width_eta, eta_range)
eta_bins_background, dsigma_background_eta = calculate_dsigma(eta_leptons_background, background_cross_section, bin_width_eta, eta_range)

pt_jet_bins_signal_0, dsigma_signal_jet_pt_0 = calculate_dsigma(pt_leading_jet_signal_0, signal_cross_section_0, bin_width_pt_jet, pt_range_jet)
pt_jet_bins_signal_2, dsigma_signal_jet_pt_2 = calculate_dsigma(pt_leading_jet_signal_2, signal_cross_section_2, bin_width_pt_jet, pt_range_jet)
pt_jet_bins_background, dsigma_background_jet_pt = calculate_dsigma(pt_leading_jet_background, background_cross_section, bin_width_pt_jet, pt_range_jet)



# Plot the differential cross-sections
# plt.figure(figsize=(10, 8))

plt.step(pt_bins_signal_0, dsigma_signal_pt_0, where="mid", alpha=0.7, label="LHeC@1.2 TeV : Signal ($w^+ w^-) [f_{M_0} / \Lambda^4$]", color="red", linewidth=3)
plt.step(pt_bins_signal_2, dsigma_signal_pt_2, where="mid", alpha=0.7, label="LHeC@1.2 TeV : Signal ($w^+ w^-) [f_{M_2} / \Lambda^4$]", color="green", linewidth=3)
plt.step(pt_bins_background, dsigma_background_pt, where="mid", alpha=0.7, label="LHeC@1.2 TeV : SM background ($w^+ w^-$)", color="blue", linewidth=3)
plt.xlabel(r"$p_T^{\ell} \ \mathrm{[GeV]}$")
plt.ylabel(r"$\frac{d\sigma}{dp_T^{\ell}} \ \mathrm{[pb/GeV]}$")
plt.title(r"$e^- p \to e^- w^+ w^- p \to e^- j j \ell \nu_{\ell} p$ : LHeC@1.2 TeV", fontsize=24)
plt.yscale("log")
plt.legend()
plt.grid(True, linestyle="--", alpha=0.6)
plt.tight_layout()
plt.ylim(0.000001, 10.0)
plt.savefig("differential_cross_section_pt.png", dpi=600)
plt.show()



# Plot the differential cross-sections for eta
#plt.figure(figsize=(10, 8))

plt.step(eta_bins_signal_0, dsigma_signal_eta_0, where="mid", alpha=0.7, label="LHeC@1.2 TeV : Signal ($w^+ w^-) [f_{M_0} / \Lambda^4$]", color="red", linewidth=3)
plt.step(eta_bins_signal_2, dsigma_signal_eta_2, where="mid", alpha=0.7, label="LHeC@1.2 TeV : Signal ($w^+ w^-) [f_{M_2} / \Lambda^4$]", color="green", linewidth=3)
plt.step(eta_bins_background, dsigma_background_eta, where="mid", alpha=0.7, label="LHeC@1.2 TeV : SM background ($w^+ w^-$)", color="blue", linewidth=3)
plt.xlabel(r"$\eta^{\ell}$")
plt.ylabel(r"$\frac{d\sigma}{d\eta^{\ell}} \ \mathrm{[pb]}$")
plt.title(r"$e^- p \to e^- w^+ w^- p \to e^- j j \ell \nu_{\ell} p$ : LHeC@1.2 TeV", fontsize=24)
plt.yscale("log")
plt.ylim(0.000001, 10.0)
plt.legend()
plt.grid(True, linestyle="--", alpha=0.6)
plt.tight_layout()
plt.savefig("differential_cross_section_eta.png", dpi=600)
plt.show()



# Plot the differential cross-sections for leading jet pT
#plt.figure(figsize=(10, 8))

plt.step(pt_jet_bins_signal_0, dsigma_signal_jet_pt_0, where="mid", alpha=0.7, label="LHeC@1.2 TeV : Signal ($w^+ w^-) [f_{M_0} / \Lambda^4$]", color="red", linewidth=3)
plt.step(pt_jet_bins_signal_2, dsigma_signal_jet_pt_2, where="mid", alpha=0.7, label="LHeC@1.2 TeV : Signal ($w^+ w^-) [f_{M_2} / \Lambda^4$]", color="green", linewidth=3)
plt.step(pt_jet_bins_background, dsigma_background_jet_pt, where="mid", alpha=0.7, label="LHeC@1.2 TeV : SM background ($w^+ w^-$)", color="blue", linewidth=3)
plt.xlabel(r"$p_T^{\mathrm{leading~jet}} \ \mathrm{[GeV]}$")
plt.ylabel(r"$\frac{d\sigma}{dp_T^{\mathrm{leading~jet}}} \ \mathrm{[pb/GeV]}$")
plt.title(r"$e^- p \to e^- w^+ w^- p \to e^- j j \ell \nu_{\ell} p$ : LHeC@1.2 TeV", fontsize=24)
plt.yscale("log")
plt.ylim(0.000001, 10.0)
plt.legend()
plt.grid(True, linestyle="--", alpha=0.6)
plt.tight_layout()
plt.savefig("differential_cross_section_leading_jet_pt.png", dpi=600)
plt.show()



# Save data for pT lepton
np.savetxt("dsigma_signal_pt.txt", np.column_stack([pt_bins_signal_0, dsigma_signal_pt_0]), header="pT [GeV], dσ/dpT [pb/GeV]")
np.savetxt("dsigma_background_pt.txt", np.column_stack([pt_bins_background, dsigma_background_pt]), header="pT [GeV], dσ/dpT [pb/GeV]")


# Save data for eta lepton
np.savetxt("dsigma_signal_eta.txt", np.column_stack([eta_bins_signal_0, dsigma_signal_eta_0]), header="eta, dσ/deta [pb]")
np.savetxt("dsigma_background_eta.txt", np.column_stack([eta_bins_background, dsigma_background_eta]), header="eta, dσ/deta [pb]")


# Save data for leading jet pT
np.savetxt("dsigma_signal_leading_jet_pt.txt", np.column_stack([pt_jet_bins_signal_0, dsigma_signal_jet_pt_0]), header="pT [GeV], dσ/dpT [pb/GeV]")
np.savetxt("dsigma_background_leading_jet_pt.txt", np.column_stack([pt_jet_bins_background, dsigma_background_jet_pt]), header="pT [GeV], dσ/dpT [pb/GeV]")


