
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



# Function to parse the LHE file and extract lepton transverse momentum (pT), leading jet pT, Delta R, Missing Transverse Energy, and Centrality

def parse_lhe_file(file_name):
    pt_leptons = []
    eta_leptons = []
    pt_leading_jet = []  # To store the transverse momentum of the leading jet
    delta_r_values = []  # To store Delta R between leptons and leading jet
    missing_transverse_energy = []  # To store Missing Transverse Energy (MET)
    centrality_values = []  # To store the centrality of the lepton
    exp_centrality_values = []  # To store the exponential centrality of the lepton
    jet_centrality_values = []  # To store the centrality of the jets
    delta_eta_jj_values = []  # To store the pseudorapidity difference between jets

    with open(file_name, "r") as file:
        in_event = False
        total_px = 0.0  # To calculate MET
        total_py = 0.0  # To calculate MET

        for line in file:
            line = line.strip()

            # Check if we are inside an <event> block
            if "<event>" in line:
                in_event = True
                jets = []  # Temporary list to store jet 4-momenta in an event
                leptons = []  # Temporary list to store lepton 4-momenta
                total_px = 0.0  # Reset for each event
                total_py = 0.0  # Reset for each event
                continue
            if "</event>" in line:
                in_event = False

                # Calculate MET (Missing Transverse Energy)
                met = np.sqrt(total_px**2 + total_py**2)
                missing_transverse_energy.append(met)

                # Determine the leading jet in this event (if any)
                if jets:
                    leading_jet = max(jets, key=lambda jet: jet.Pt())  # Find the leading jet
                    pt_leading_jet.append(leading_jet.Pt())  # Store the transverse momentum of the leading jet

                    # Calculate Delta R between each lepton and the leading jet
                    for lepton in leptons:
                        delta_r = lepton.DeltaR(leading_jet)
                        delta_r_values.append(delta_r)

                    # Calculate lepton centrality if there are at least two jets
                    if len(jets) >= 2 and len(leptons) == 1:
                        delta_eta_jj = abs(jets[0].Eta() - jets[1].Eta())
                        delta_eta_jj_values.append(delta_eta_jj)  # Store Δηjj
                        if delta_eta_jj > 0:  # Avoid division by zero
                            lepton_centrality = (leptons[0].Eta() - jets[0].Eta()) / delta_eta_jj
                            centrality_values.append(lepton_centrality)

                    # Calculate the centrality of jets
                    if len(jets) >= 2:  # Ensure there are at least two jets
                        jet_centrality = abs(jets[0].Eta() + jets[1].Eta()) / 2.0
                        jet_centrality_values.append(jet_centrality)

                # Calculate the new centrality (e^(-eta) / cosh(eta)) for each lepton
                for lepton in leptons:
                    exp_centrality = np.exp(-lepton.Eta()) / np.cosh(lepton.Eta())
                    exp_centrality_values.append(exp_centrality)

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

                # Accumulate px and py for MET calculation
                total_px += px
                total_py += py

                # Check if the particle is a lepton (e- or e+)
                if abs(pdg_id) in [11, 13, 15]:  # e, mu, tau
                    # Create TLorentzVector
                    lepton = ROOT.TLorentzVector()
                    lepton.SetPxPyPzE(px, py, pz, energy)

                    # Store lepton 4-momentum
                    leptons.append(lepton)

                    # Extract transverse momentum and pseudorapidity
                    pt_leptons.append(lepton.Pt())
                    eta_leptons.append(lepton.Eta())

                # Check if the particle is a jet (status 1, generic criterion for outgoing particles)
                if abs(pdg_id) not in [11, 12, 13, 14, 15, 16, 22] and int(parts[1]) == 1:
                    # Create TLorentzVector for jet
                    jet = ROOT.TLorentzVector()
                    jet.SetPxPyPzE(px, py, pz, energy)

                    # Store jet 4-momentum
                    jets.append(jet)

    return (pt_leptons, eta_leptons, pt_leading_jet, delta_r_values,
            missing_transverse_energy, centrality_values, exp_centrality_values,
            jet_centrality_values, delta_eta_jj_values)



#fig, ax = plt.subplots(figsize=(12.0, 10.0))
plt.subplots_adjust(left=0.15, right=0.95, bottom=0.12, top=0.95)



# Parse signal and background files
signal_file_0 = "/home/hamzeh-khanpour/MG5_aMC_v3_5_7/aa_ww_semi_leptonic_NP_FM0/Events/run_01/aa_ww_semi_leptonic_NP_FM0.lhe"
signal_file_2 = "/home/hamzeh-khanpour/MG5_aMC_v3_5_7/aa_ww_semi_leptonic_NP_FM2/Events/run_01/aa_ww_semi_leptonic_NP_FM2.lhe"
background_file = "/home/hamzeh-khanpour/MG5_aMC_v3_5_7/aa_ww_semi_leptonic_SM/Events/run_01/aa_ww_semi_leptonic_SM.lhe"





# Parse signal and background files, including Delta R, MET, centrality, exp_centrality, jet_centrality, and Delta_eta_jj
(pt_leptons_signal_0, eta_leptons_signal_0, pt_leading_jet_signal_0, delta_r_signal_0,
 met_signal_0, centrality_signal_0, exp_centrality_signal_0, jet_centrality_signal_0, delta_eta_jj_signal_0) = parse_lhe_file(signal_file_0)

(pt_leptons_signal_2, eta_leptons_signal_2, pt_leading_jet_signal_2, delta_r_signal_2,
 met_signal_2, centrality_signal_2, exp_centrality_signal_2, jet_centrality_signal_2, delta_eta_jj_signal_2) = parse_lhe_file(signal_file_2)

(pt_leptons_background, eta_leptons_background, pt_leading_jet_background, delta_r_background,
 met_background, centrality_background, exp_centrality_background, jet_centrality_background, delta_eta_jj_background) = parse_lhe_file(background_file)





# Parameters for differential cross-section
signal_cross_section_0 = 1.994  # pb
signal_cross_section_2 = 77.86  # pb
background_cross_section = 0.0134802  # pb
num_bins = 50
pt_range_lepton = (0, 500)     # Range for lepton pT
pt_range_jet = (0, 500)        # Range for leading jet pT (adjusted for higher jet momenta)
eta_range = (-10, 10)          # Range for pseudorapidity
delta_r_range = (0, 10)        # Range for Delta R is between 0 and 5
met_range = (1, 500)           # Define range for MET (adjust as needed)
centrality_range = (-5, 5)     # Centrality typically ranges from -5 to 5
exp_centrality_range = (0, 2)  # Centrality typically ranges from 0 to 2
jet_centrality_range = (0, 6) # Centrality typically ranges of jet from 0 to 10
delta_eta_jj_range  = (0, 5)   # Pseudorapidity difference between jets from 0 to 5




# Calculate bin width
bin_width_pt_lepton = (pt_range_lepton[1] - pt_range_lepton[0]) / num_bins
bin_width_pt_jet = (pt_range_jet[1] - pt_range_jet[0]) / num_bins
bin_width_eta = (eta_range[1] - eta_range[0]) / num_bins
bin_width_delta_r = (delta_r_range[1] - delta_r_range[0]) / num_bins
bin_width_met = (met_range[1] - met_range[0]) / num_bins  # Bin width for MET
bin_width_centrality = (centrality_range[1] - centrality_range[0]) / num_bins  # Bin width for lepton centrality
bin_width_exp_centrality = (exp_centrality_range[1] - exp_centrality_range[0]) / num_bins  # Bin width for exponential centrality
bin_width_jet_centrality = (jet_centrality_range[1] - jet_centrality_range[0]) / num_bins  # Bin width for jet centrality
bin_width_delta_eta_jj = (delta_eta_jj_range[1] - delta_eta_jj_range[0]) / num_bins  # Bin width for Δηjj





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


# Normalize Delta R histograms
delta_r_bins_signal_0, dsigma_signal_delta_r_0 = calculate_dsigma(delta_r_signal_0, signal_cross_section_0, bin_width_delta_r, delta_r_range)
delta_r_bins_signal_2, dsigma_signal_delta_r_2 = calculate_dsigma(delta_r_signal_2, signal_cross_section_2, bin_width_delta_r, delta_r_range)
delta_r_bins_background, dsigma_background_delta_r = calculate_dsigma(delta_r_background, background_cross_section, bin_width_delta_r, delta_r_range)


# Calculate differential cross-sections for MET
met_bins_signal_0, dsigma_signal_met_0 = calculate_dsigma(met_signal_0, signal_cross_section_0, bin_width_met, met_range)
met_bins_signal_2, dsigma_signal_met_2 = calculate_dsigma(met_signal_2, signal_cross_section_2, bin_width_met, met_range)
met_bins_background, dsigma_background_met = calculate_dsigma(met_background, background_cross_section, bin_width_met, met_range)



# Calculate differential cross-sections for Lepton Centrality
centrality_bins_signal_0, dsigma_signal_centrality_0 = calculate_dsigma(centrality_signal_0, signal_cross_section_0, bin_width_centrality, centrality_range)
centrality_bins_signal_2, dsigma_signal_centrality_2 = calculate_dsigma(centrality_signal_2, signal_cross_section_2, bin_width_centrality, centrality_range)
centrality_bins_background, dsigma_background_centrality = calculate_dsigma(centrality_background, background_cross_section, bin_width_centrality, centrality_range)



# Calculate differential cross-sections for exp_centrality
exp_centrality_bins_signal_0, dsigma_signal_exp_centrality_0 = calculate_dsigma(exp_centrality_signal_0, signal_cross_section_0, bin_width_exp_centrality, centrality_range)
exp_centrality_bins_signal_2, dsigma_signal_exp_centrality_2 = calculate_dsigma(exp_centrality_signal_2, signal_cross_section_2, bin_width_exp_centrality, centrality_range)
exp_centrality_bins_background, dsigma_background_exp_centrality = calculate_dsigma(exp_centrality_background, background_cross_section, bin_width_exp_centrality, centrality_range)



# Calculate differential cross-sections for jet_centrality
jet_centrality_bins_signal_0, dsigma_signal_jet_centrality_0 = calculate_dsigma(jet_centrality_signal_0, signal_cross_section_0, bin_width_jet_centrality, jet_centrality_range)
jet_centrality_bins_signal_2, dsigma_signal_jet_centrality_2 = calculate_dsigma(jet_centrality_signal_2, signal_cross_section_2, bin_width_jet_centrality, jet_centrality_range)
jet_centrality_bins_background, dsigma_background_jet_centrality = calculate_dsigma(jet_centrality_background, background_cross_section, bin_width_jet_centrality, jet_centrality_range)



# Calculate differential cross-sections for Δηjj
delta_eta_jj_bins_signal_0, dsigma_signal_delta_eta_jj_0 = calculate_dsigma(delta_eta_jj_signal_0, signal_cross_section_0, bin_width_delta_eta_jj, delta_eta_jj_range)
delta_eta_jj_bins_signal_2, dsigma_signal_delta_eta_jj_2 = calculate_dsigma(delta_eta_jj_signal_2, signal_cross_section_2, bin_width_delta_eta_jj, delta_eta_jj_range)
delta_eta_jj_bins_background, dsigma_background_delta_eta_jj = calculate_dsigma(delta_eta_jj_background, background_cross_section, bin_width_delta_eta_jj, delta_eta_jj_range)








# =======================================================================


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
plt.ylim(0.00001, 10.0)
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
plt.ylim(0.0001, 100.0)
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
plt.ylim(0.00001, 10.0)
plt.legend()
plt.grid(True, linestyle="--", alpha=0.6)
plt.tight_layout()
plt.savefig("differential_cross_section_leading_jet_pt.png", dpi=600)
plt.show()




# Plot Delta R differential cross-section
plt.step(delta_r_bins_signal_0, dsigma_signal_delta_r_0, where="mid", alpha=0.7, label="LHeC@1.2 TeV : Signal ($w^+ w^-) [f_{M_0} / \Lambda^4$]", color="red", linewidth=3)
plt.step(delta_r_bins_signal_2, dsigma_signal_delta_r_2, where="mid", alpha=0.7, label="LHeC@1.2 TeV : Signal ($w^+ w^-) [f_{M_2} / \Lambda^4$]", color="green", linewidth=3)
plt.step(delta_r_bins_background, dsigma_background_delta_r, where="mid", alpha=0.7, label="LHeC@1.2 TeV : SM background ($w^+ w^-$)", color="blue", linewidth=3)
plt.xlabel(r"$\Delta R(\ell, \mathrm{leading~jet})$")
plt.ylabel(r"$\frac{d\sigma}{d\Delta R} \ \mathrm{[pb]}$")
plt.title(r"$e^- p \to e^- w^+ w^- p \to e^- j j \ell \nu_{\ell} p$ : LHeC@1.2 TeV", fontsize=24)
plt.yscale("log")
plt.ylim(0.0001, 1000.0)
plt.legend()
plt.grid(True, linestyle="--", alpha=0.6)
plt.tight_layout()
plt.savefig("differential_cross_section_delta_r.png", dpi=600)
plt.show()





# Plot the differential cross-sections for Missing Transverse Energy (MET)
plt.step(met_bins_signal_0, dsigma_signal_met_0, where="mid", alpha=0.7, label="LHeC@1.2 TeV : Signal ($w^+ w^-) [f_{M_0} / \Lambda^4$]", color="red", linewidth=3)
plt.step(met_bins_signal_2, dsigma_signal_met_2, where="mid", alpha=0.7, label="LHeC@1.2 TeV : Signal ($w^+ w^-) [f_{M_2} / \Lambda^4$]", color="green", linewidth=3)
plt.step(met_bins_background, dsigma_background_met, where="mid", alpha=0.7, label="LHeC@1.2 TeV : SM background ($w^+ w^-$)", color="blue", linewidth=3)
plt.xlabel(r"$\mathrm{MET} \ \mathrm{[GeV]}$")
plt.ylabel(r"$\frac{d\sigma}{d\mathrm{MET}} \ \mathrm{[pb/GeV]}$")
plt.title(r"$e^- p \to e^- w^+ w^- p \to e^- j j \ell \nu_{\ell} p$ : LHeC@1.2 TeV", fontsize=24)
plt.yscale("log")
plt.legend()
plt.grid(True, linestyle="--", alpha=0.6)
plt.tight_layout()
plt.ylim(0.000001, 1.0)
plt.savefig("differential_cross_section_met.png", dpi=600)
plt.show()







# Normalize the distributions
dsigma_signal_centrality_0_norm = dsigma_signal_centrality_0 / np.sum(dsigma_signal_centrality_0)
dsigma_signal_centrality_2_norm = dsigma_signal_centrality_2 / np.sum(dsigma_signal_centrality_2)
dsigma_background_centrality_norm = dsigma_background_centrality / np.sum(dsigma_background_centrality)

# Plot the normalized distributions for Lepton Centrality
plt.step(centrality_bins_signal_0, dsigma_signal_centrality_0_norm, where="mid", alpha=0.7, label="LHeC@1.2 TeV : Signal ($w^+ w^-) [f_{M_0} / \Lambda^4$]", color="red", linewidth=3)
plt.step(centrality_bins_signal_2, dsigma_signal_centrality_2_norm, where="mid", alpha=0.7, label="LHeC@1.2 TeV : Signal ($w^+ w^-) [f_{M_2} / \Lambda^4$]", color="green", linewidth=3)
plt.step(centrality_bins_background, dsigma_background_centrality_norm, where="mid", alpha=0.7, label="LHeC@1.2 TeV : SM background ($w^+ w^-$)", color="blue", linewidth=3)
# Axis labels and title
plt.xlabel(r"$C_{\ell}$")
plt.ylabel("Normalized Distribution")
plt.title(r"$e^- p \to e^- w^+ w^- p \to e^- j j \ell \nu_{\ell} p$ : LHeC@1.2 TeV", fontsize=24)
#plt.yscale("log")
plt.ylim(0.0, 0.08)  # Adjust as needed for log scale
# Add legend, grid, and formula
plt.legend()
plt.grid(True, linestyle="--", alpha=0.6)
plt.text(0.5, 1e-3, r"$C_{\ell} = \frac{\eta_{\ell} - \frac{\eta_{\mathrm{jet1}} + \eta_{\mathrm{jet2}}}{2}}{\Delta \eta_{jj}}$", color="black")
# Save and display the plot
plt.tight_layout()
plt.savefig("normalized_distribution_centrality.png", dpi=600)
plt.show()



# Plot the differential cross-sections for Lepton Centrality
plt.step(centrality_bins_signal_0, dsigma_signal_centrality_0, where="mid", alpha=0.7, label="LHeC@1.2 TeV : Signal ($w^+ w^-) [f_{M_0} / \Lambda^4$]", color="red", linewidth=3)
plt.step(centrality_bins_signal_2, dsigma_signal_centrality_2, where="mid", alpha=0.7, label="LHeC@1.2 TeV : Signal ($w^+ w^-) [f_{M_2} / \Lambda^4$]", color="green", linewidth=3)
plt.step(centrality_bins_background, dsigma_background_centrality, where="mid", alpha=0.7, label="LHeC@1.2 TeV : SM background ($w^+ w^-$)", color="blue", linewidth=3)
plt.xlabel(r"$C_{\ell}$")
plt.ylabel(r"$\frac{d\sigma}{dC_{\ell}} \ \mathrm{[pb]}$")
plt.title(r"$e^- p \to e^- w^+ w^- p \to e^- j j \ell \nu_{\ell} p$ : LHeC@1.2 TeV", fontsize=24)
plt.yscale("log")
plt.ylim(0.00001, 100.0)
plt.legend()
plt.grid(True, linestyle="--", alpha=0.6)
# Add formula inside the plot
plt.text(0.5, 0.001, r"$C_{\ell} = \frac{\eta_{\ell} - \frac{\eta_{\mathrm{jet1}} + \eta_{\mathrm{jet2}}}{2}}{\Delta \eta_{jj}}$",  color="black")
plt.tight_layout()
plt.savefig("differential_cross_section_centrality.png", dpi=600)
plt.show()









# Normalize the distributions
dsigma_signal_exp_centrality_0_norm = dsigma_signal_exp_centrality_0 / np.sum(dsigma_signal_exp_centrality_0)
dsigma_signal_exp_centrality_2_norm = dsigma_signal_exp_centrality_2 / np.sum(dsigma_signal_exp_centrality_2)
dsigma_background_exp_centrality_norm = dsigma_background_exp_centrality / np.sum(dsigma_background_exp_centrality)

# Plot the normalized distributions for Exponential Lepton Centrality
plt.step(exp_centrality_bins_signal_0, dsigma_signal_exp_centrality_0_norm, where="mid", alpha=0.7, label="LHeC@1.2 TeV : Signal ($w^+ w^-) [f_{M_0} / \Lambda^4$)", color="red", linewidth=3)
plt.step(exp_centrality_bins_signal_2, dsigma_signal_exp_centrality_2_norm, where="mid", alpha=0.7, label="LHeC@1.2 TeV : Signal ($w^+ w^-) [f_{M_2} / \Lambda^4$)", color="green", linewidth=3)
plt.step(exp_centrality_bins_background, dsigma_background_exp_centrality_norm, where="mid", alpha=0.7, label="LHeC@1.2 TeV : SM background ($w^+ w^-$)", color="blue", linewidth=3)
# Axis labels and title
plt.xlabel(r"$C_{\ell}^{\mathrm{exp}}$")
plt.ylabel("Normalized Distribution")
plt.title(r"$e^- p \to e^- w^+ w^- p \to e^- j j \ell \nu_{\ell} p$ : LHeC@1.2 TeV", fontsize=24)
plt.yscale("log")
plt.ylim(1e-5, 1.0)  # Adjust as needed for log scale
# Add legend, grid, and formula
plt.legend()
plt.grid(True, linestyle="--", alpha=0.6)
plt.text(0.5, 1e-3, r"$C_{\ell}^{\mathrm{exp}} = e^{-|C_{\ell}|}$", color="black")
# Save and display the plot
plt.tight_layout()
plt.savefig("normalized_distribution_exp_centrality.png", dpi=600)
plt.show()




# Plot the differential cross-sections for Exponential Lepton Centrality (exp_centrality)
plt.step(exp_centrality_bins_signal_0, dsigma_signal_exp_centrality_0, where="mid", alpha=0.7, label="LHeC@1.2 TeV : Signal ($w^+ w^-) [f_{M_0} / \Lambda^4$]", color="red", linewidth=3)
plt.step(exp_centrality_bins_signal_2, dsigma_signal_exp_centrality_2, where="mid", alpha=0.7, label="LHeC@1.2 TeV : Signal ($w^+ w^-) [f_{M_2} / \Lambda^4$]", color="green", linewidth=3)
plt.step(exp_centrality_bins_background, dsigma_background_exp_centrality, where="mid", alpha=0.7, label="LHeC@1.2 TeV : SM background ($w^+ w^-$)", color="blue", linewidth=3)
plt.xlabel(r"$C_{\ell}^{\mathrm{exp}}$")
plt.ylabel(r"$\frac{d\sigma}{dC_{\ell}^{\mathrm{exp}}} \ \mathrm{[pb]}$")
plt.title(r"$e^- p \to e^- w^+ w^- p \to e^- j j \ell \nu_{\ell} p$ : LHeC@1.2 TeV", fontsize=24)
plt.yscale("log")
plt.ylim(0.00001, 100.0)
plt.legend()
plt.grid(True, linestyle="--", alpha=0.6)
# Add formula inside the plot
plt.text(0.5, 0.001, r"$C_{\ell}^{\mathrm{exp}} = e^{-|C_{\ell}|}$", color="black")
plt.tight_layout()
plt.savefig("differential_cross_section_exp_centrality.png", dpi=600)
plt.show()











# Normalize the distributions
dsigma_signal_jet_centrality_0_norm = dsigma_signal_jet_centrality_0 / np.sum(dsigma_signal_jet_centrality_0)
dsigma_signal_jet_centrality_2_norm = dsigma_signal_jet_centrality_2 / np.sum(dsigma_signal_jet_centrality_2)
dsigma_background_jet_centrality_norm = dsigma_background_jet_centrality / np.sum(dsigma_background_jet_centrality)

# Plot the normalized distributions for Jet Centrality
plt.step(jet_centrality_bins_signal_0, dsigma_signal_jet_centrality_0_norm, where="mid", alpha=0.7, label="LHeC@1.2 TeV : Signal ($w^+ w^-) [f_{M_0} / \Lambda^4$)", color="red", linewidth=3)
plt.step(jet_centrality_bins_signal_2, dsigma_signal_jet_centrality_2_norm, where="mid", alpha=0.7, label="LHeC@1.2 TeV : Signal ($w^+ w^-) [f_{M_2} / \Lambda^4$)", color="green", linewidth=3)
plt.step(jet_centrality_bins_background, dsigma_background_jet_centrality_norm, where="mid", alpha=0.7, label="LHeC@1.2 TeV : SM background ($w^+ w^-$)", color="blue", linewidth=3)
# Axis labels and title
plt.xlabel(r"$C_{\mathrm{jets}}$")
plt.ylabel("Normalized Distribution")
plt.title(r"$e^- p \to e^- w^+ w^- p \to e^- j j \ell \nu_{\ell} p$ : LHeC@1.2 TeV", fontsize=24)
#plt.yscale("log")
plt.ylim(0.0, 0.08)  # Adjusted for normalized scale
# Add legend, grid, and formula inside the plot
plt.legend()
plt.grid(True, linestyle="--", alpha=0.6)
plt.text(0.5, 1e-3, r"$C_{\mathrm{jets}} = \frac{|\eta_{\mathrm{jet1}} + \eta_{\mathrm{jet2}}|}{2}$", color="black")
# Save and display the plot
plt.tight_layout()
plt.savefig("normalized_distribution_jet_centrality.png", dpi=600)
plt.show()



# Plot the differential cross-sections for Jet Centrality
plt.step(jet_centrality_bins_signal_0, dsigma_signal_jet_centrality_0, where="mid", alpha=0.7, label="LHeC@1.2 TeV : Signal ($w^+ w^-) [f_{M_0} / \Lambda^4$)", color="red", linewidth=3)
plt.step(jet_centrality_bins_signal_2, dsigma_signal_jet_centrality_2, where="mid", alpha=0.7, label="LHeC@1.2 TeV : Signal ($w^+ w^-) [f_{M_2} / \Lambda^4$)", color="green", linewidth=3)
plt.step(jet_centrality_bins_background, dsigma_background_jet_centrality, where="mid", alpha=0.7, label="LHeC@1.2 TeV : SM background ($w^+ w^-$)", color="blue", linewidth=3)
plt.xlabel(r"$C_{\mathrm{jets}}$")
plt.ylabel(r"$\frac{d\sigma}{dC_{\mathrm{jets}}} \ \mathrm{[pb]}$")
plt.title(r"$e^- p \to e^- w^+ w^- p \to e^- j j \ell \nu_{\ell} p$ : LHeC@1.2 TeV", fontsize=24)
plt.yscale("log")
plt.ylim(0.00001, 100.0)
plt.legend()
plt.grid(True, linestyle="--", alpha=0.6)
# Add formula inside the plot
plt.text(0.5, 0.001, r"$C_{\mathrm{jets}} = \frac{|\eta_{\mathrm{jet1}} + \eta_{\mathrm{jet2}}|}{2}$", color="black")
plt.tight_layout()
plt.savefig("differential_cross_section_jet_centrality.png", dpi=600)
plt.show()













# Plot the differential cross-sections for Pseudorapidity Difference Between Jets (Δηjj)
plt.step(delta_eta_jj_bins_signal_0, dsigma_signal_delta_eta_jj_0, where="mid", alpha=0.7, label="LHeC@1.2 TeV : Signal ($w^+ w^-) [f_{M_0} / \Lambda^4$)", color="red", linewidth=3)
plt.step(delta_eta_jj_bins_signal_2, dsigma_signal_delta_eta_jj_2, where="mid", alpha=0.7, label="LHeC@1.2 TeV : Signal ($w^+ w^-) [f_{M_2} / \Lambda^4$)", color="green", linewidth=3)
plt.step(delta_eta_jj_bins_background, dsigma_background_delta_eta_jj, where="mid", alpha=0.7, label="LHeC@1.2 TeV : SM background ($w^+ w^-$)", color="blue", linewidth=3)
plt.xlabel(r"$\Delta \eta_{jj}$")
plt.ylabel(r"$\frac{d\sigma}{d\Delta \eta_{jj}} \ \mathrm{[pb]}$")
plt.title(r"$e^- p \to e^- w^+ w^- p \to e^- j j \ell \nu_{\ell} p$ : LHeC@1.2 TeV", fontsize=24)
plt.yscale("log")
plt.ylim(0.0001, 1000.0)
plt.legend()
plt.grid(True, linestyle="--", alpha=0.6)
plt.tight_layout()
plt.savefig("differential_cross_section_delta_eta_jj.png", dpi=600)
plt.show()







# =======================================================================






# Save data for pT lepton
np.savetxt("dsigma_signal_pt_0.txt", np.column_stack([pt_bins_signal_0, dsigma_signal_pt_0]), header="pT [GeV], dσ/dpT [pb/GeV]")
np.savetxt("dsigma_signal_pt_2.txt", np.column_stack([pt_bins_signal_2, dsigma_signal_pt_2]), header="pT [GeV], dσ/dpT [pb/GeV]")
np.savetxt("dsigma_background_pt.txt", np.column_stack([pt_bins_background, dsigma_background_pt]), header="pT [GeV], dσ/dpT [pb/GeV]")


# Save data for eta lepton
np.savetxt("dsigma_signal_eta_0.txt", np.column_stack([eta_bins_signal_0, dsigma_signal_eta_0]), header="eta, dσ/deta [pb]")
np.savetxt("dsigma_signal_eta_2.txt", np.column_stack([eta_bins_signal_2, dsigma_signal_eta_2]), header="eta, dσ/deta [pb]")
np.savetxt("dsigma_background_eta.txt", np.column_stack([eta_bins_background, dsigma_background_eta]), header="eta, dσ/deta [pb]")


# Save data for leading jet pT
np.savetxt("dsigma_signal_leading_jet_pt_0.txt", np.column_stack([pt_jet_bins_signal_0, dsigma_signal_jet_pt_0]), header="pT [GeV], dσ/dpT [pb/GeV]")
np.savetxt("dsigma_signal_leading_jet_pt_2.txt", np.column_stack([pt_jet_bins_signal_2, dsigma_signal_jet_pt_2]), header="pT [GeV], dσ/dpT [pb/GeV]")
np.savetxt("dsigma_background_leading_jet_pt.txt", np.column_stack([pt_jet_bins_background, dsigma_background_jet_pt]), header="pT [GeV], dσ/dpT [pb/GeV]")


# Save data for Delta R
np.savetxt("dsigma_signal_delta_r_0.txt", np.column_stack([delta_r_bins_signal_0, dsigma_signal_delta_r_0]), header="Delta R, dσ/d(Delta R) [pb]")
np.savetxt("dsigma_signal_delta_r_2.txt", np.column_stack([delta_r_bins_signal_2, dsigma_signal_delta_r_2]), header="Delta R, dσ/d(Delta R) [pb]")
np.savetxt("dsigma_background_delta_r.txt", np.column_stack([delta_r_bins_background, dsigma_background_delta_r]), header="Delta R, dσ/d(Delta R) [pb]")


# Save data for Missing Transverse Energy (MET)
np.savetxt("dsigma_signal_met_0.txt", np.column_stack([met_bins_signal_0, dsigma_signal_met_0]), header="MET [GeV], dσ/d(MET) [pb/GeV]")
np.savetxt("dsigma_signal_met_2.txt", np.column_stack([met_bins_signal_2, dsigma_signal_met_2]), header="MET [GeV], dσ/d(MET) [pb/GeV]")
np.savetxt("dsigma_background_met.txt", np.column_stack([met_bins_background, dsigma_background_met]), header="MET [GeV], dσ/d(MET) [pb/GeV]")


# Save data for Lepton Centrality
np.savetxt("dsigma_signal_centrality_0.txt", np.column_stack([centrality_bins_signal_0, dsigma_signal_centrality_0]), header="C_l, dσ/d(C_l) [pb]")
np.savetxt("dsigma_signal_centrality_2.txt", np.column_stack([centrality_bins_signal_2, dsigma_signal_centrality_2]), header="C_l, dσ/d(C_l) [pb]")
np.savetxt("dsigma_background_centrality.txt", np.column_stack([centrality_bins_background, dsigma_background_centrality]), header="C_l, dσ/d(C_l) [pb]")


# Save data for Exponential Lepton Centrality (exp_centrality)
np.savetxt("dsigma_signal_exp_centrality_0.txt", np.column_stack([exp_centrality_bins_signal_0, dsigma_signal_exp_centrality_0]), header="C_l_exp, dσ/d(C_l_exp) [pb]")
np.savetxt("dsigma_signal_exp_centrality_2.txt", np.column_stack([exp_centrality_bins_signal_2, dsigma_signal_exp_centrality_2]), header="C_l_exp, dσ/d(C_l_exp) [pb]")
np.savetxt("dsigma_background_exp_centrality.txt", np.column_stack([exp_centrality_bins_background, dsigma_background_exp_centrality]), header="C_l_exp, dσ/d(C_l_exp) [pb]")


# Save data for Jet Centrality
np.savetxt("dsigma_signal_jet_centrality_0.txt", np.column_stack([jet_centrality_bins_signal_0, dsigma_signal_jet_centrality_0]), header="C_j, dσ/d(C_j) [pb]")
np.savetxt("dsigma_signal_jet_centrality_2.txt", np.column_stack([jet_centrality_bins_signal_2, dsigma_signal_jet_centrality_2]), header="C_j, dσ/d(C_j) [pb]")
np.savetxt("dsigma_background_jet_centrality.txt", np.column_stack([jet_centrality_bins_background, dsigma_background_jet_centrality]), header="C_j, dσ/d(C_j) [pb]")



# Save data for Pseudorapidity Difference Between Jets (Δηjj)
np.savetxt("dsigma_signal_delta_eta_jj_0.txt", np.column_stack([delta_eta_jj_bins_signal_0, dsigma_signal_delta_eta_jj_0]), header="Δη_jj, dσ/d(Δη_jj) [pb]")
np.savetxt("dsigma_signal_delta_eta_jj_2.txt", np.column_stack([delta_eta_jj_bins_signal_2, dsigma_signal_delta_eta_jj_2]), header="Δη_jj, dσ/d(Δη_jj) [pb]")
np.savetxt("dsigma_background_delta_eta_jj.txt", np.column_stack([delta_eta_jj_bins_background, dsigma_background_delta_eta_jj]), header="Δη_jj, dσ/d(Δη_jj) [pb]")







