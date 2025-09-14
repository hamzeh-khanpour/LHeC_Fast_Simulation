
# =======================================================================
#     Hamzeh khanpour -- 2025   -- LHeC
# =======================================================================


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


# =======================================================================


# Function to parse the LHE file and extract lepton transverse momentum (pT), leading jet pT, Delta R, Missing Transverse Energy, and Centrality
def parse_lhe_file(file_name):
    pt_leptons = []
    eta_leptons = []
    pt_leading_jet = []  # Leading jet pT
    delta_r_values = []  # Delta R between leptons and leading jet
    missing_transverse_energy = []  # MET (Missing Transverse Energy)
    centrality_values = []  # Lepton centrality
    exp_centrality_values = []  # Exponential centrality
    jet_centrality_values = []  # Jet centrality
    delta_eta_jj_values = []  # Pseudorapidity difference between jets
    m_w_hadronic_values = []  # Hadronic W boson mass (W → jj)
    m_w_leptonic_values = []  # Leptonic W boson mass (W → lν)

    m_w_hadronic_leptonic_values = []  # WW invariant mass


    with open(file_name, "r") as file:
        in_event = False

        for line in file:
            line = line.strip()

            # Check if we are inside an <event> block
            if "<event>" in line:
                in_event = True
                jets = []  # Temporary list to store jet 4-momenta
                leptons = []  # Temporary list to store lepton 4-momenta
                neutrinos = []  # Store neutrino 4-momenta for MET calculation
                continue
            if "</event>" in line:
                in_event = False

                # Compute MET from neutrinos' pT
                met_x = sum(nu.Px() for nu in neutrinos)
                met_y = sum(nu.Py() for nu in neutrinos)
                met = np.sqrt(met_x**2 + met_y**2)
                missing_transverse_energy.append(met)

                # Identify the leading jet in this event (if available)
                if jets:
                    leading_jet = max(jets, key=lambda jet: jet.Pt())  # Select highest pT jet
                    pt_leading_jet.append(leading_jet.Pt())  # Store leading jet pT

                    # Compute ΔR between each lepton and the leading jet
                    for lepton in leptons:
                        delta_r = lepton.DeltaR(leading_jet)
                        delta_r_values.append(delta_r)

                    # Compute lepton centrality if at least two jets exist
                    if len(jets) >= 2 and len(leptons) == 1:
                        delta_eta_jj = jets[0].Eta()
                        delta_eta_jj_values.append(delta_eta_jj)  # Store Δηjj
                        if delta_eta_jj > -10000:  # Avoid division by zero
                            lepton_centrality = jets[1].Eta()
                            centrality_values.append(lepton_centrality)

                    # Compute the centrality of jets (average pseudorapidity of two jets)
                    if len(jets) >= 2:
                        jet_centrality = abs(jets[0].Eta() + jets[1].Eta()) / 2.0
                        jet_centrality_values.append(jet_centrality)

                    # Compute hadronic W boson mass if at least two jets exist
                    if len(jets) >= 2:
                        w_hadronic = jets[0] + jets[1]  # Sum of two jets (4-vectors)
                        m_w_hadronic_values.append(w_hadronic.M())  # Store mass

                # Compute leptonic W boson mass if one lepton and at least one neutrino exist
                if len(leptons) == 1 and len(neutrinos) > 0:
                    neutrino_vec = sum(neutrinos, ROOT.TLorentzVector())  # Sum of all neutrinos
                    w_leptonic = leptons[0] + neutrino_vec
                    m_w_leptonic_values.append(w_leptonic.M())  # Store mass

                # --- M_WW (reco) ---
                if w_hadronic is not None and w_leptonic is not None:
                    w_hadronic_leptonic = w_hadronic + w_leptonic
                    m_w_hadronic_leptonic_values.append(w_hadronic_leptonic.M())

                # Compute the exponential centrality for each lepton
                for lepton in leptons:
                    exp_centrality = np.exp(-lepton.Eta()) / np.cosh(lepton.Eta())
                    exp_centrality_values.append(exp_centrality)

                continue

            if in_event:
                # Ignore metadata lines (only parse numerical data)
                if not line[0].isdigit() and not line[0] == "-":
                    continue

                # Extract numerical data
                parts = line.split()
                if len(parts) < 10:
                    continue

                pdg_id = int(parts[0])
                px = float(parts[6])
                py = float(parts[7])
                pz = float(parts[8])
                energy = float(parts[9])

                # Check if particle is a lepton (electron, muon, tau)
                if abs(pdg_id) in [11, 13, 15]:
                    lepton = ROOT.TLorentzVector()
                    lepton.SetPxPyPzE(px, py, pz, energy)

                    # Store lepton momentum and properties
                    leptons.append(lepton)
                    pt_leptons.append(lepton.Pt())
                    eta_leptons.append(lepton.Eta())

                # Check if the particle is a jet (not a lepton, neutrino, or photon)
                if abs(pdg_id) not in [11, 12, 13, 14, 15, 16, 22] and int(parts[1]) == 1:
                    jet = ROOT.TLorentzVector()
                    jet.SetPxPyPzE(px, py, pz, energy)

                    # Store jet information
                    jets.append(jet)

                # Check if the particle is a neutrino (νe, νμ, ντ)
                if abs(pdg_id) in [12, 14, 16]:  # Neutrinos
                    neutrino = ROOT.TLorentzVector()
                    neutrino.SetPxPyPzE(px, py, pz, energy)

                    # Store neutrino 4-momentum for MET calculation
                    neutrinos.append(neutrino)



    return (pt_leptons, eta_leptons, pt_leading_jet, delta_r_values,
            missing_transverse_energy, centrality_values, exp_centrality_values,
            jet_centrality_values, delta_eta_jj_values, m_w_hadronic_values, m_w_leptonic_values, m_w_hadronic_leptonic_values)






# =======================================================================





#fig, ax = plt.subplots(figsize=(12.0, 10.0))
plt.subplots_adjust(left=0.15, right=0.95, bottom=0.12, top=0.95)



# Parse signal and background files
signal_file_0 = "/home/hamzeh-khanpour/MG5_aMC_v3_6_3/aa_ww_semi_leptonic_NP_1_FM2/Events/run_02/aa_ww_semi_leptonic_NP_1_FM2_10.lhe"
#signal_file_1 = "/home/hamzeh-khanpour/MG5_aMC_v3_5_7/aa_ww_semi_leptonic_NP_1_FM1/Events/run_01/aa_ww_semi_leptonic_NP_1_FM1.lhe"
signal_file_2 = "/home/hamzeh-khanpour/MG5_aMC_v3_6_3/aa_ww_semi_leptonic_NP_1_FM2/Events/run_03/aa_ww_semi_leptonic_NP_1_FM2_-10.lhe"
#signal_file_3 = "/home/hamzeh-khanpour/MG5_aMC_v3_5_7/aa_ww_semi_leptonic_NP_1-FM3/Events/run_01/aa_ww_semi_leptonic_NP_1_FM3.lhe"


background_file = "/home/hamzeh-khanpour/MG5_aMC_v3_6_3/aa_ww_semi_leptonic_SM_NP_0_FMi_0/Events/run_01/aa_ww_semi_leptonic_SM_NP_0_FMi_0.lhe"




# Parse signal and background files, including Delta R, MET, centrality, exp_centrality, jet_centrality, and W masses (m_w_hadronic, m_w_leptonic)
# Signal FM0
(pt_leptons_signal_0, eta_leptons_signal_0, pt_leading_jet_signal_0, delta_r_signal_0,
 met_signal_0, centrality_signal_0, exp_centrality_signal_0, jet_centrality_signal_0, delta_eta_jj_signal_0,
 m_w_hadronic_signal_0, m_w_leptonic_signal_0, m_w_hadronic_leptonic_signal_0) = parse_lhe_file(signal_file_0)

# Signal FM2
(pt_leptons_signal_2, eta_leptons_signal_2, pt_leading_jet_signal_2, delta_r_signal_2,
 met_signal_2, centrality_signal_2, exp_centrality_signal_2, jet_centrality_signal_2, delta_eta_jj_signal_2,
 m_w_hadronic_signal_2, m_w_leptonic_signal_2, m_w_hadronic_leptonic_signal_2) = parse_lhe_file(signal_file_2)

# SM background
(pt_leptons_background, eta_leptons_background, pt_leading_jet_background, delta_r_background,
 met_background, centrality_background, exp_centrality_background, jet_centrality_background, delta_eta_jj_background,
 m_w_hadronic_background, m_w_leptonic_background, m_w_hadronic_leptonic_background) = parse_lhe_file(background_file)






# =======================================================================



# Parameters for differential cross-section




# 100  NP = 1
#signal_cross_section_0    = 0.01490319      # pb  FM0
#signal_cross_section_1   = 0.01508150      # pb  FM1
#signal_cross_section_2    = 0.021425        # pb  FM2
#signal_cross_section_3   = 0.01644609      # pb  FM3



# 10   NP = 1
signal_cross_section_0    = 0.0153508       # pb  FM2 = 10
#signal_cross_section_1   = 0.0             # pb  FM1
signal_cross_section_2    = 0.015646        # pb  FM2 = -10
#signal_cross_section_3   = 0.0             # pb  FM3


background_cross_section  = 0.0149219    # pb     real : 0.0154074




num_bins = 40

pt_range_lepton = (0, 250)     # Range for lepton pT
pt_range_jet = (0, 250)        # Range for leading jet pT (adjusted for higher jet momenta)
eta_range = (-10, 10)          # Range for pseudorapidity
delta_r_range = (0, 10)        # Range for Delta R is between 0 and 5
met_range = (1, 250)           # Define range for MET (adjust as needed)
centrality_range = (-10, 10)     # Centrality typically ranges from -5 to 5
exp_centrality_range = (0, 2)  # Centrality typically ranges from 0 to 2
jet_centrality_range = (0, 6) # Centrality typically ranges of jet from 0 to 10
delta_eta_jj_range  = (-10, 10)   # Pseudorapidity difference between jets from 0 to 5
m_w_hadronic_range = (1, 140)  # Range for the hadronic W boson mass
m_w_leptonic_range = (1, 140)  # Range for the leptonic W boson mass

# >>> ADD THESE LINES FOR RECO M_WW <<<
m_w_hadronic_leptonic_range = (160, 1000)  # WW invariant mass; threshold ~ 2*mW up to ~LHeC √s_γγ





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
# Calculate bin widths for various distributions
bin_width_m_w_hadronic = (m_w_hadronic_range[1] - m_w_hadronic_range[0]) / num_bins  # Bin width for M_W hadronic
bin_width_m_w_leptonic = (m_w_leptonic_range[1] - m_w_leptonic_range[0]) / num_bins  # Bin width for M_W leptonic

bin_width_m_w_hadronic_leptonic = (m_w_hadronic_leptonic_range[1] - m_w_hadronic_leptonic_range[0]) / num_bins







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



# Calculate differential cross-sections for Hadronic W boson mass (M_W^jj)
m_w_hadronic_bins_signal_0, dsigma_signal_m_w_hadronic_0 = calculate_dsigma(m_w_hadronic_signal_0, signal_cross_section_0, bin_width_m_w_hadronic, m_w_hadronic_range)
m_w_hadronic_bins_signal_2, dsigma_signal_m_w_hadronic_2 = calculate_dsigma(m_w_hadronic_signal_2, signal_cross_section_2, bin_width_m_w_hadronic, m_w_hadronic_range)
m_w_hadronic_bins_background, dsigma_background_m_w_hadronic = calculate_dsigma(m_w_hadronic_background, background_cross_section, bin_width_m_w_hadronic, m_w_hadronic_range)



# Calculate differential cross-sections for Leptonic W boson mass (M_W^ℓν)
m_w_leptonic_bins_signal_0, dsigma_signal_m_w_leptonic_0 = calculate_dsigma(m_w_leptonic_signal_0, signal_cross_section_0, bin_width_m_w_leptonic, m_w_leptonic_range)
m_w_leptonic_bins_signal_2, dsigma_signal_m_w_leptonic_2 = calculate_dsigma(m_w_leptonic_signal_2, signal_cross_section_2, bin_width_m_w_leptonic, m_w_leptonic_range)
m_w_leptonic_bins_background, dsigma_background_m_w_leptonic = calculate_dsigma(m_w_leptonic_background, background_cross_section, bin_width_m_w_leptonic, m_w_leptonic_range)



# >>> ADD THESE LINES FOR RECO M_WW <<<
m_w_hadronic_leptonic_bins_signal_0, dsigma_signal_m_w_hadronic_leptonic_0 = calculate_dsigma(    m_w_hadronic_leptonic_signal_0, signal_cross_section_0, bin_width_m_w_hadronic_leptonic, m_w_hadronic_leptonic_range)
m_w_hadronic_leptonic_bins_signal_2, dsigma_signal_m_w_hadronic_leptonic_2 = calculate_dsigma(    m_w_hadronic_leptonic_signal_2, signal_cross_section_2, bin_width_m_w_hadronic_leptonic, m_w_hadronic_leptonic_range)
m_w_hadronic_leptonic_bins_background, dsigma_background_m_w_hadronic_leptonic = calculate_dsigma(    m_w_hadronic_leptonic_background, background_cross_section, bin_width_m_w_hadronic_leptonic, m_w_hadronic_leptonic_range)








# =======================================================================




# Plot the differential cross-sections
# plt.figure(figsize=(10, 8))

plt.step(pt_bins_signal_0, dsigma_signal_pt_0, where="mid", alpha=0.7, label="LHeC@1.2 TeV : Signal ($w^+ w^-) [f_{M_1} / \Lambda^4=10$]", color="red", linewidth=3)
plt.step(pt_bins_signal_2, dsigma_signal_pt_2, where="mid", alpha=0.7, label="LHeC@1.2 TeV : Signal ($w^+ w^-) [f_{M_2} / \Lambda^4=-10$]", color="green", linewidth=3)
plt.step(pt_bins_background, dsigma_background_pt, where="mid", alpha=0.7, label="LHeC@1.2 TeV : SM background ($w^+ w^-$)", color="blue", linewidth=3)
plt.xlabel(r"$p_T^{\ell} \ \mathrm{[GeV]}$")
plt.ylabel(r"$\frac{d\sigma}{dp_T^{\ell}} \ \mathrm{[pb/GeV]}$")
plt.title(r"$e^- p \to e^- w^+ w^- p \to e^- j j \ell \nu_{\ell} p$ : LHeC@1.2 TeV", fontsize=24)
plt.yscale("log")
plt.legend()
plt.grid(True, linestyle="--", alpha=0.6)
plt.tight_layout()
plt.ylim(0.000001, 0.01)
plt.savefig("differential_cross_section_pt_LHeC.pdf", dpi=600)
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
plt.ylim(0.00001, 0.1)
plt.legend()
plt.grid(True, linestyle="--", alpha=0.6)
plt.tight_layout()
#plt.savefig("differential_cross_section_eta_lepton_LHeC.pdf", dpi=600)
plt.show()




# Plot the differential cross-sections for leading jet pT
#plt.figure(figsize=(10, 8))

plt.step(pt_jet_bins_signal_0, dsigma_signal_jet_pt_0, where="mid", alpha=0.7, label="LHeC@1.2 TeV : Signal ($w^+ w^-) [f_{M_0} / \Lambda^4=10$]", color="red", linewidth=3)
plt.step(pt_jet_bins_signal_2, dsigma_signal_jet_pt_2, where="mid", alpha=0.7, label="LHeC@1.2 TeV : Signal ($w^+ w^-) [f_{M_2} / \Lambda^4=-10$]", color="green", linewidth=3)
plt.step(pt_jet_bins_background, dsigma_background_jet_pt, where="mid", alpha=0.7, label="LHeC@1.2 TeV : SM background ($w^+ w^-$)", color="blue", linewidth=3)
plt.xlabel(r"$p_T^{\mathrm{leading~jet}} \ \mathrm{[GeV]}$")
plt.ylabel(r"$\frac{d\sigma}{dp_T^{\mathrm{leading~jet}}} \ \mathrm{[pb/GeV]}$")
plt.title(r"$e^- p \to e^- w^+ w^- p \to e^- j j \ell \nu_{\ell} p$ : LHeC@1.2 TeV", fontsize=24)
plt.yscale("log")
plt.ylim(0.000001, 0.01)
plt.legend()
plt.grid(True, linestyle="--", alpha=0.6)
plt.tight_layout()
plt.savefig("differential_cross_section_leading_jet_pt_LHeC.pdf", dpi=600)
plt.show()






# Plot Delta R differential cross-section
plt.step(delta_r_bins_signal_0, dsigma_signal_delta_r_0, where="mid", alpha=0.7, label="LHeC@1.2 TeV : Signal ($w^+ w^-) [f_{M_0} / \Lambda^4$]", color="red", linewidth=3)
plt.step(delta_r_bins_signal_2, dsigma_signal_delta_r_2, where="mid", alpha=0.7, label="LHeC@1.2 TeV : Signal ($w^+ w^-) [f_{M_2} / \Lambda^4$]", color="green", linewidth=3)
plt.step(delta_r_bins_background, dsigma_background_delta_r, where="mid", alpha=0.7, label="LHeC@1.2 TeV : SM background ($w^+ w^-$)", color="blue", linewidth=3)
plt.xlabel(r"$\Delta R(\ell, \mathrm{leading~jet})$")
plt.ylabel(r"$\frac{d\sigma}{d\Delta R} \ \mathrm{[pb]}$")
plt.title(r"$e^- p \to e^- w^+ w^- p \to e^- j j \ell \nu_{\ell} p$ : LHeC@1.2 TeV", fontsize=24)
plt.yscale("log")
plt.ylim(0.0001, 1.0)
plt.legend()
plt.grid(True, linestyle="--", alpha=0.6)
plt.tight_layout()
#plt.savefig("differential_cross_section_delta_r_LHeC.pdf", dpi=600)
plt.show()





# Plot the differential cross-sections for Missing Transverse Energy (MET)
plt.step(met_bins_signal_0, dsigma_signal_met_0, where="mid", alpha=0.7, label="LHeC@1.2 TeV : Signal ($w^+ w^-) [f_{M_0} / \Lambda^4=10$]", color="red", linewidth=3)
plt.step(met_bins_signal_2, dsigma_signal_met_2, where="mid", alpha=0.7, label="LHeC@1.2 TeV : Signal ($w^+ w^-) [f_{M_2} / \Lambda^4=-10$]", color="green", linewidth=3)
plt.step(met_bins_background, dsigma_background_met, where="mid", alpha=0.7, label="LHeC@1.2 TeV : SM background ($w^+ w^-$)", color="blue", linewidth=3)
plt.xlabel(r"$\mathrm{MET} \ \mathrm{[GeV]}$")
plt.ylabel(r"$\frac{d\sigma}{d\mathrm{MET}} \ \mathrm{[pb/GeV]}$")
plt.title(r"$e^- p \to e^- w^+ w^- p \to e^- j j \ell \nu_{\ell} p$ : LHeC@1.2 TeV", fontsize=24)
plt.yscale("log")
plt.legend()
plt.grid(True, linestyle="--", alpha=0.6)
plt.tight_layout()
plt.ylim(0.000001, 0.01)
plt.savefig("differential_cross_section_met_LHeC.pdf", dpi=600)
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
plt.xlabel(r"$\eta_{leading \; jet}$")
plt.ylabel("Normalized Distribution")
plt.title(r"$e^- p \to e^- w^+ w^- p \to e^- j j \ell \nu_{\ell} p$ : LHeC@1.2 TeV", fontsize=24)
#plt.yscale("log")
plt.ylim(0.0, 0.2)  # Adjust as needed for log scale
# Add legend, grid, and formula
plt.legend()
plt.grid(True, linestyle="--", alpha=0.6)
#plt.text(0.5, 1e-3, r"$C_{\ell} = \frac{\eta_{\ell} - \frac{\eta_{\mathrm{jet1}} + \eta_{\mathrm{jet2}}}{2}}{\Delta \eta_{jj}}$", color="black")
# Save and display the plot
plt.tight_layout()
#plt.savefig("normalized_distribution_eta_leading_jet_LHeC.pdf", dpi=600)
plt.show()






# Plot the differential cross-sections for Lepton Centrality
plt.step(centrality_bins_signal_0, dsigma_signal_centrality_0, where="mid", alpha=0.7, label="LHeC@1.2 TeV : Signal ($w^+ w^-) [f_{M_0} / \Lambda^4$]", color="red", linewidth=3)
plt.step(centrality_bins_signal_2, dsigma_signal_centrality_2, where="mid", alpha=0.7, label="LHeC@1.2 TeV : Signal ($w^+ w^-) [f_{M_2} / \Lambda^4$]", color="green", linewidth=3)
plt.step(centrality_bins_background, dsigma_background_centrality, where="mid", alpha=0.7, label="LHeC@1.2 TeV : SM background ($w^+ w^-$)", color="blue", linewidth=3)
plt.xlabel(r"$\eta_{leading \; jet}$")
plt.ylabel(r"$\frac{d\sigma}{d\eta_{leading \; jet}}} \ \mathrm{[pb]}$")
plt.title(r"$e^- p \to e^- w^+ w^- p \to e^- j j \ell \nu_{\ell} p$ : LHeC@1.2 TeV", fontsize=24)
plt.yscale("log")
plt.ylim(0.00001, 0.1)
plt.legend()
plt.grid(True, linestyle="--", alpha=0.6)
# Add formula inside the plot
#plt.text(0.5, 0.001, r"$C_{\ell} = \frac{\eta_{\ell} - \frac{\eta_{\mathrm{jet1}} + \eta_{\mathrm{jet2}}}{2}}{\Delta \eta_{jj}}$",  color="black")
plt.tight_layout()
#plt.savefig("differential_cross_section_eta_leading_jet_LHeC.pdf", dpi=600)
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
#plt.savefig("normalized_distribution_exp_centrality_LHeC.pdf", dpi=600)
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
#plt.savefig("differential_cross_section_exp_centrality_LHeC.pdf", dpi=600)
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
plt.ylim(0.02, 0.08)  # Adjusted for normalized scale
# Add legend, grid, and formula inside the plot
plt.legend()
plt.grid(True, linestyle="--", alpha=0.6)
plt.text(0.9, 0.06, r"$C_{\mathrm{jets}} = \frac{|\eta_{\mathrm{jet1}} + \eta_{\mathrm{jet2}}|}{2}$", color="black")
# Save and display the plot
plt.tight_layout()
#plt.savefig("normalized_distribution_jet_centrality_LHeC.pdf", dpi=600)
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
#plt.savefig("differential_cross_section_jet_centrality_LHeC.pdf", dpi=600)
plt.show()








# Plot the differential cross-sections for Pseudorapidity Difference Between Jets (Δηjj)
plt.step(delta_eta_jj_bins_signal_0, dsigma_signal_delta_eta_jj_0, where="mid", alpha=0.7, label="LHeC@1.2 TeV : Signal ($w^+ w^-) [f_{M_0} / \Lambda^4$)", color="red", linewidth=3)
plt.step(delta_eta_jj_bins_signal_2, dsigma_signal_delta_eta_jj_2, where="mid", alpha=0.7, label="LHeC@1.2 TeV : Signal ($w^+ w^-) [f_{M_2} / \Lambda^4$)", color="green", linewidth=3)
plt.step(delta_eta_jj_bins_background, dsigma_background_delta_eta_jj, where="mid", alpha=0.7, label="LHeC@1.2 TeV : SM background ($w^+ w^-$)", color="blue", linewidth=3)
plt.xlabel(r"$\eta_{subleading \; jet}$")
plt.ylabel(r"$\frac{d\sigma}{d\eta_{subleading \; jet}} \ \mathrm{[pb]}$")
plt.title(r"$e^- p \to e^- w^+ w^- p \to e^- j j \ell \nu_{\ell} p$ : LHeC@1.2 TeV", fontsize=24)
plt.yscale("log")
plt.ylim(0.0001, 0.1)
plt.legend()
plt.grid(True, linestyle="--", alpha=0.6)
plt.tight_layout()
#plt.savefig("differential_cross_section_eta_subleading_jet_LHeC.pdf", dpi=600)
plt.show()





# Normalize the hadronic W boson mass distributions
dsigma_signal_m_w_hadronic_0_norm = dsigma_signal_m_w_hadronic_0 / np.sum(dsigma_signal_m_w_hadronic_0)
dsigma_signal_m_w_hadronic_2_norm = dsigma_signal_m_w_hadronic_2 / np.sum(dsigma_signal_m_w_hadronic_2)
dsigma_background_m_w_hadronic_norm = dsigma_background_m_w_hadronic / np.sum(dsigma_background_m_w_hadronic)



# --------------------------------------------
# Plot the normalized Hadronic W Boson Mass (M_W^jj)
#plt.figure(figsize=(10, 8))

plt.step(m_w_hadronic_bins_signal_0, dsigma_signal_m_w_hadronic_0_norm, where="mid", alpha=0.7, label="LHeC@1.2 TeV : Signal ($w^+ w^-) [f_{M_0} / \Lambda^4$)", color="red", linewidth=3)
plt.step(m_w_hadronic_bins_signal_2, dsigma_signal_m_w_hadronic_2_norm, where="mid", alpha=0.7, label="LHeC@1.2 TeV : Signal ($w^+ w^-) [f_{M_2} / \Lambda^4$)", color="green", linewidth=3)
plt.step(m_w_hadronic_bins_background, dsigma_background_m_w_hadronic_norm, where="mid", alpha=0.7, label="LHeC@1.2 TeV : SM background ($w^+ w^-$)", color="blue", linewidth=3)
plt.xlabel(r"$M_W^{jj} \ \mathrm{[GeV]}$")
plt.ylabel(r"Normalized Distribution")
plt.title(r"$e^- p \to e^- w^+ w^- p \to e^- j j \ell \nu_{\ell} p$ : LHeC@1.2 TeV", fontsize=24)
#plt.yscale("log")
#plt.ylim(0.0, 0.08)  # Adjusted for normalized scale
plt.legend()
plt.grid(True, linestyle="--", alpha=0.6)
plt.tight_layout()
#plt.savefig("normalized_cross_section_m_w_hadronic_LHeC.pdf", dpi=600)
plt.show()






# --------------------------------------------
# Plot the differential cross-sections for Hadronic W Boson Mass (M_W^jj)
plt.step(m_w_hadronic_bins_signal_0, dsigma_signal_m_w_hadronic_0, where="mid", alpha=0.7, label="LHeC@1.2 TeV : Signal ($w^+ w^-) [f_{M_0} / \Lambda^4$)", color="red", linewidth=3)
plt.step(m_w_hadronic_bins_signal_2, dsigma_signal_m_w_hadronic_2, where="mid", alpha=0.7, label="LHeC@1.2 TeV : Signal ($w^+ w^-) [f_{M_2} / \Lambda^4$)", color="green", linewidth=3)
plt.step(m_w_hadronic_bins_background, dsigma_background_m_w_hadronic, where="mid", alpha=0.7, label="LHeC@1.2 TeV : SM background ($w^+ w^-$)", color="blue", linewidth=3)
plt.xlabel(r"$M_W^{jj} \ \mathrm{[GeV]}$")
plt.ylabel(r"$\frac{d\sigma}{dM_W^{jj}} \ \mathrm{[pb/GeV]}$")
plt.title(r"$e^- p \to e^- w^+ w^- p \to e^- j j \ell \nu_{\ell} p$ : LHeC@1.2 TeV", fontsize=24)
plt.yscale("log")
plt.ylim(0.00001, 0.1)
plt.legend()
plt.grid(True, linestyle="--", alpha=0.6)
plt.tight_layout()
#plt.savefig("differential_cross_section_m_w_hadronic_LHeC.pdf", dpi=600)
plt.show()





# Normalize the leptonic W boson mass distributions
dsigma_signal_m_w_leptonic_0_norm = dsigma_signal_m_w_leptonic_0 / np.sum(dsigma_signal_m_w_leptonic_0)
dsigma_signal_m_w_leptonic_2_norm = dsigma_signal_m_w_leptonic_2 / np.sum(dsigma_signal_m_w_leptonic_2)
dsigma_background_m_w_leptonic_norm = dsigma_background_m_w_leptonic / np.sum(dsigma_background_m_w_leptonic)





# --------------------------------------------
# Plot the normalized Leptonic W Boson Mass (M_W^ℓν)
#plt.figure(figsize=(10, 8))
plt.step(m_w_leptonic_bins_signal_0, dsigma_signal_m_w_leptonic_0_norm, where="mid", alpha=0.7, label="LHeC@1.2 TeV : Signal ($w^+ w^-) [f_{M_0} / \Lambda^4$)", color="red", linewidth=3)
plt.step(m_w_leptonic_bins_signal_2, dsigma_signal_m_w_leptonic_2_norm, where="mid", alpha=0.7, label="LHeC@1.2 TeV : Signal ($w^+ w^-) [f_{M_2} / \Lambda^4$)", color="green", linewidth=3)
plt.step(m_w_leptonic_bins_background, dsigma_background_m_w_leptonic_norm, where="mid", alpha=0.7, label="LHeC@1.2 TeV : SM background ($w^+ w^-$)", color="blue", linewidth=3)
plt.xlabel(r"$M_W^{\ell\nu_{\ell}} \ \mathrm{[GeV]}$")
plt.ylabel(r"Normalized Distribution")
plt.title(r"$e^- p \to e^- w^+ w^- p \to e^- j j \ell \nu_{\ell} p$ : LHeC@1.2 TeV", fontsize=24)
#plt.yscale("log")
#plt.ylim(0.0, 0.08)  # Adjusted for normalized scale
plt.legend()
plt.grid(True, linestyle="--", alpha=0.6)
plt.tight_layout()
#plt.savefig("normalized_cross_section_m_w_leptonic_LHeC.pdf", dpi=600)
plt.show()






# --------------------------------------------
# Plot the differential cross-sections for Leptonic W Boson Mass (M_W^ℓν)
plt.step(m_w_leptonic_bins_signal_0, dsigma_signal_m_w_leptonic_0, where="mid", alpha=0.7, label="LHeC@1.2 TeV : Signal ($w^+ w^-) [f_{M_0} / \Lambda^4$)", color="red", linewidth=3)
plt.step(m_w_leptonic_bins_signal_2, dsigma_signal_m_w_leptonic_2, where="mid", alpha=0.7, label="LHeC@1.2 TeV : Signal ($w^+ w^-) [f_{M_2} / \Lambda^4$)", color="green", linewidth=3)
plt.step(m_w_leptonic_bins_background, dsigma_background_m_w_leptonic, where="mid", alpha=0.7, label="LHeC@1.2 TeV : SM background ($w^+ w^-$)", color="blue", linewidth=3)
plt.xlabel(r"$M_W^{\ell\nu_{\ell}} \ \mathrm{[GeV]}$")
plt.ylabel(r"$\frac{d\sigma}{dM_W^{\ell\nu_{\ell}}} \ \mathrm{[pb/GeV]}$")
plt.title(r"$e^- p \to e^- w^+ w^- p \to e^- j j \ell \nu_{\ell} p$ : LHeC@1.2 TeV", fontsize=24)
plt.yscale("log")
plt.ylim(0.00001, 0.1)
plt.legend()
plt.grid(True, linestyle="--", alpha=0.6)
plt.tight_layout()
#plt.savefig("differential_cross_section_m_w_leptonic_LHeC.pdf", dpi=600)
plt.show()







# --- Normalize the reconstructed M_WW distributions ---
s0 = np.sum(dsigma_signal_m_w_hadronic_leptonic_0)
s2 = np.sum(dsigma_signal_m_w_hadronic_leptonic_2)
sb = np.sum(dsigma_background_m_w_hadronic_leptonic)

dsigma_signal_m_w_hadronic_leptonic_0_norm = dsigma_signal_m_w_hadronic_leptonic_0 / s0 if s0 > 0 else dsigma_signal_m_w_hadronic_leptonic_0
dsigma_signal_m_w_hadronic_leptonic_2_norm = dsigma_signal_m_w_hadronic_leptonic_2 / s2 if s2 > 0 else dsigma_signal_m_w_hadronic_leptonic_2
dsigma_background_m_w_hadronic_leptonic_norm = dsigma_background_m_w_hadronic_leptonic / sb if sb > 0 else dsigma_background_m_w_hadronic_leptonic




# --------------------------------------------
# Plot the normalized Reconstructed WW Invariant Mass (M_WW)
plt.step(m_w_hadronic_leptonic_bins_signal_0, dsigma_signal_m_w_hadronic_leptonic_0_norm, where="mid", alpha=0.7, label="LHeC@1.2 TeV : Signal ($w^+ w^-) [f_{M_2} / \\Lambda^4=10$)", color="red", linewidth=3)
plt.step(m_w_hadronic_leptonic_bins_signal_2, dsigma_signal_m_w_hadronic_leptonic_2_norm, where="mid", alpha=0.7, label="LHeC@1.2 TeV : Signal ($w^+ w^-) [f_{M_2} / \\Lambda^4=-10$)", color="green", linewidth=3)
plt.step(m_w_hadronic_leptonic_bins_background, dsigma_background_m_w_hadronic_leptonic_norm, where="mid", alpha=0.7, label="LHeC@1.2 TeV : SM background ($w^+ w^-$)", color="blue", linewidth=3)
plt.xlabel(r"$M_{WW} \ \mathrm{[GeV]}$")
plt.ylabel("Normalized Distribution")
plt.title(r"$e^- p \to e^- w^+ w^- p \to e^- j j \ell \nu_{\ell} p$ : LHeC@1.2 TeV", fontsize=24)
# plt.yscale("log")  # usually off for normalized shapes; keep commented like your leptonic plot
plt.legend()
plt.grid(True, linestyle="--", alpha=0.6)
plt.tight_layout()
plt.savefig("normalized_cross_section_m_ww_reco_LHeC.pdf", dpi=600)
plt.show()





# --------------------------------------------
# Plot the differential cross-sections for Reconstructed WW Invariant Mass (M_WW)
plt.step(m_w_hadronic_leptonic_bins_signal_0, dsigma_signal_m_w_hadronic_leptonic_0, where="mid", alpha=0.7, label="LHeC@1.2 TeV : Signal ($w^+ w^-) [f_{M_2} / \\Lambda^4=10$)", color="red", linewidth=3)
plt.step(m_w_hadronic_leptonic_bins_signal_2, dsigma_signal_m_w_hadronic_leptonic_2, where="mid", alpha=0.7, label="LHeC@1.2 TeV : Signal ($w^+ w^-) [f_{M_2} / \\Lambda^4=-10$)", color="green", linewidth=3)
plt.step(m_w_hadronic_leptonic_bins_background, dsigma_background_m_w_hadronic_leptonic, where="mid", alpha=0.7, label="LHeC@1.2 TeV : SM background ($w^+ w^-$)", color="blue", linewidth=3)
plt.xlabel(r"$M_{WW} \ \mathrm{[GeV]}$")
plt.ylabel(r"$\frac{d\sigma}{dM_{WW}} \ \mathrm{[pb/GeV]}$")
plt.title(r"$e^- p \to e^- w^+ w^- p \to e^- j j \ell \nu_{\ell} p$ : LHeC@1.2 TeV", fontsize=24)
plt.yscale("log")
plt.ylim(0.000001, 0.01)
plt.legend()
plt.grid(True, linestyle="--", alpha=0.6)
plt.tight_layout()
plt.savefig("differential_cross_section_m_ww_reco_LHeC.pdf", dpi=600)
plt.show()



#=====================================================================================



# =============================================================================
# Tail σ(MWW>W_cut) and 95% CL limits (counting, Asimov exact) using f=0, +10, −10
# =============================================================================

import math
import numpy as np


# --- Config
THRESHOLDS = [400.0, 500.0, 600.0, 700.0, 800.0]  # GeV
L_fb = 1000.0     # integrated luminosity in fb^-1
Aeps = 1.0        # acceptance × efficiency
PB_TO_FB = 1e3    # 1 pb = 1000 fb
f0 = 10.0         # TeV^-4  (your MG couplings are ±10)



# --- Build the bin edges for MWW from your (range, nbins)
edges_mww = np.linspace(m_w_hadronic_leptonic_range[0], m_w_hadronic_leptonic_range[1], num_bins + 1)


def integrate_tail_hist(y, edges, thr):
    """Integrate a density histogram y (pb/GeV) above threshold thr using bin edges."""
    total = 0.0
    for i in range(len(y)):
        lo, hi = edges[i], edges[i+1]
        if hi <= thr:
            continue
        width = hi - lo
        if lo < thr < hi:
            total += y[i] * (hi - thr)
        else:
            total += y[i] * width
    return float(total)  # pb


def q_asimov(mu, mu0):
    """One-bin Asimov test statistic (no nuisances)."""
    if mu <= 0 or mu0 <= 0:
        return float("inf")
    return 2.0 * (mu - mu0 + mu0 * math.log(mu0 / mu))


def interval_asimov_exact(c, A, B, L_fb, Aeps, q_target=3.84, f_scan=300.0):
    """
    Solve q(f)=q_target for two-sided 95% CL via bracketing+bisection.
    σ(f)=c + A f + B f^2  (pb);  μ(f) = L_fb * Aeps * σ(f) * 1000.
    """
    def qf(f):
        sig = max(c + A*f + B*f*f, 1e-18)
        mu  = L_fb * Aeps * sig * PB_TO_FB
        mu0 = L_fb * Aeps * c   * PB_TO_FB
        return q_asimov(mu, mu0)

    # Right root (f >= 0)
    flo, fhi, step = 0.0, 0.0, 0.1
    while True:
        fhi = flo + step
        if qf(fhi) >= q_target or fhi >= f_scan:
            break
        step *= 1.8
    for _ in range(80):
        mid = 0.5*(flo + fhi)
        if qf(mid) < q_target: flo = mid
        else:                  fhi = mid
    f_plus = 0.5*(flo + fhi)

    # Left root (f <= 0)
    fhi, flo, step = 0.0, 0.0, -0.1
    while True:
        flo = fhi + step
        if qf(flo) >= q_target or abs(flo) >= f_scan:
            break
        step *= 1.8
    for _ in range(80):
        mid = 0.5*(flo + fhi)
        if qf(mid) < q_target: fhi = mid
        else:                  flo = mid
    f_minus = 0.5*(flo + fhi)

    return f_minus, f_plus


def extract_AB(c, sp, sm, f0):
    """
    From tail cross sections at f=0, +f0, −f0 (pb), build:
      σ(f) = c + A f + B f^2
    """
    A = (sp - sm) / (2.0 * f0)
    B = (sp + sm - 2.0 * c) / (2.0 * f0 * f0)
    return A, B



# --- Tail σ rows for each cut (pb)
row_sm = [integrate_tail_hist(dsigma_background_m_w_hadronic_leptonic, edges_mww, thr) for thr in THRESHOLDS]  # f=0
row_p  = [integrate_tail_hist(dsigma_signal_m_w_hadronic_leptonic_0, edges_mww, thr)   for thr in THRESHOLDS]  # f=+10
row_m  = [integrate_tail_hist(dsigma_signal_m_w_hadronic_leptonic_2, edges_mww, thr)   for thr in THRESHOLDS]  # f=−10


# --- Pretty print the tail σ and expected counts
print("\n=== σ(MWW > W_cut) from LHE-level spectra (pb) ===")
hdr = "cut  ".ljust(8) + "σ_SM      σ(+10)    σ(−10)    ".ljust(32) + "N_SM@1000fb  N(+10)@1000fb  N(−10)@1000fb"
print(hdr)
print("-"*len(hdr))
for thr, c, sp, sm in zip(THRESHOLDS, row_sm, row_p, row_m):
    N_c  = L_fb * PB_TO_FB * Aeps * c
    N_sp = L_fb * PB_TO_FB * Aeps * sp
    N_sm = L_fb * PB_TO_FB * Aeps * sm
    print(f"W>{int(thr):<4} {c:9.3e} {sp:9.3e} {sm:9.3e}   {N_c:12.1f}   {N_sp:12.1f}    {N_sm:12.1f}")



# --- Limits per cut; choose best (narrowest) interval
print("\n=== 95% CL (two-sided, Asimov SM) on f (TeV^-4), per cut and best ===")
print("cut   A[pb/TeV^4]       B[pb/TeV^8]        interval [f−, f+]  ")
best = None
for thr, c, sp, sm in zip(THRESHOLDS, row_sm, row_p, row_m):
    A, B = extract_AB(c, sp, sm, f0)
    f_lo, f_hi = interval_asimov_exact(c, A, B, L_fb=L_fb, Aeps=Aeps, q_target=3.84)
    width = f_hi - f_lo
    print(f"W>{int(thr):<4} {A:14.3e}  {B:14.3e}   [{f_lo:7.3g}, {f_hi:7.3g}]")
    if (best is None) or (width < best[-1]):
        best = (thr, f_lo, f_hi, A, B, c, width)

thr_best, f_lo, f_hi, A_best, B_best, c_best, _ = best
print(f"\n>> Best cut: W > {thr_best:.0f} GeV")
print(f"   σ_tail(f) = c + A f + B f^2 with c={c_best:.4e} pb, A={A_best:.4e} pb/TeV^4, B={B_best:.4e} pb/TeV^8")
print(f"   Two-sided 95% CL: f ∈ [{f_lo:.3g}, {f_hi:.3g}] TeV^-4")
# =============================================================================





# =============================================================================
# Multi-bin SHAPE likelihood in MWW (Asimov SM, no extra backgrounds)
# Uses all bins between Wmin and Wmax (full spectrum if Wmin=0.0, Wmax=None)
# =============================================================================
import numpy as np
import math

# --- Choose which ranges you want to report (full + a few tails):
SHAPE_CUTS = [(0.0, None), (400.0, None), (600.0, None), (800.0, None), (1000.0, None), (1200.0, None), (1400.0, None)]



# Build MWW bin edges from your config
edges_mww = np.linspace(m_w_hadronic_leptonic_range[0],
                        m_w_hadronic_leptonic_range[1],
                        num_bins + 1)

def _split_first_bin_if_needed(edges, y, Wmin):
    """If Wmin cuts through the first used bin, scale that bin's width proportionally."""
    lo = edges[:-1]; hi = edges[1:]
    mask = hi > Wmin
    if not np.any(mask):
        return np.array([]), np.array([])
    widths = hi - lo
    widths = widths[mask].copy()
    ysel   = np.asarray(y)[mask].copy()
    # partial coverage on the very first selected bin
    i0 = np.argmax(mask)  # index in full arrays of first selected bin
    if lo[i0] < Wmin < hi[i0]:
        frac = (hi[i0] - Wmin) / (hi[i0] - lo[i0])
        widths[0] *= frac  # keep density the same, reduce width
    return ysel, widths



def build_shape_coeffs(edges, y_sm, y_p, y_m, f0, Wmin=0.0, Wmax=None):
    """Return per-bin α,β,γ (pb/GeV) and widths (GeV) in the chosen MWW window."""
    if Wmax is None:
        Wmax = edges[-1]
    lo = edges[:-1]; hi = edges[1:]
    # pre-mask by window, then handle first-partial-bin correction
    mask_window = (hi > Wmin) & (lo < Wmax)
    edges_w = np.r_[lo[mask_window], hi[mask_window][-1]]  # local edges (not used except for widths)
    # take subarrays then correct first-bin width if partial
    α_full = np.asarray(y_sm)[mask_window]
    β_full = (np.asarray(y_p)[mask_window] - np.asarray(y_m)[mask_window]) / (2.0 * f0)
    γ_full = (np.asarray(y_p)[mask_window] + np.asarray(y_m)[mask_window] - 2.0*np.asarray(y_sm)[mask_window]) / (2.0 * f0 * f0)

    widths_full = (hi - lo)[mask_window].copy()
    # correct first bin if Wmin cuts inside
    i0 = np.argmax(mask_window)  # first used bin index in full histogram
    if lo[i0] < Wmin < hi[i0]:
        frac = (hi[i0] - Wmin) / (hi[i0] - lo[i0])
        widths_full[0] *= frac

    # clip densities to avoid negative μ from tiny MC fluctuations
    α = np.clip(α_full, 1e-18, None)
    β = β_full.copy()
    γ = γ_full.copy()
    return α, β, γ, widths_full



def asimov_shape_interval(α, β, γ, widths, L_fb, Aeps, q_target=3.84, f_scan=1e3):
    """
    Multi-bin Asimov q(f) = 2 sum_i [ μ_i(f) - μ_i(0) + μ_i(0) ln( μ_i(0) / μ_i(f) ) ]
    with μ_i(f) = L*Aε*1000 * (α_i + β_i f + γ_i f^2) * ΔW_i.
    Returns two-sided 95% CL [f-, f+].
    """
    const = L_fb * Aeps * 1000.0
    n0 = const * α * widths  # Asimov observed (SM only)

    def q_of(f):
        dens = α + β*f + γ*f*f
        dens = np.clip(dens, 1e-18, None)
        mu   = const * dens * widths
        return 2.0 * np.sum(mu - n0 + n0 * np.log(n0 / mu))

    # --- right root (f >= 0)
    flo, fhi, step = 0.0, 0.0, 0.1
    while True:
        fhi = flo + step
        if q_of(fhi) >= q_target or fhi >= f_scan:
            break
        step *= 1.8
    for _ in range(80):
        mid = 0.5*(flo+fhi)
        if q_of(mid) < q_target: flo = mid
        else: fhi = mid
    f_plus = 0.5*(flo+fhi)

    # --- left root (f <= 0)
    fhi, flo, step = 0.0, 0.0, -0.1
    while True:
        flo = fhi + step
        if q_of(flo) >= q_target or abs(flo) >= f_scan:
            break
        step *= 1.8
    for _ in range(80):
        mid = 0.5*(flo+fhi)
        if q_of(mid) < q_target: fhi = mid
        else: flo = mid
    f_minus = 0.5*(flo+fhi)
    return f_minus, f_plus



# --- Run shape limits for your requested windows
print("\n=== SHAPE likelihood limits in MWW (Asimov SM, two-sided 95% CL) ===")
best_shape = None
for (Wmin, Wmax) in SHAPE_CUTS:
    α, β, γ, widths = build_shape_coeffs(edges_mww,
                                         dsigma_background_m_w_hadronic_leptonic,
                                         dsigma_signal_m_w_hadronic_leptonic_0,  # +10
                                         dsigma_signal_m_w_hadronic_leptonic_2,  # -10
                                         f0=f0, Wmin=Wmin, Wmax=Wmax)
    if α.size == 0:
        continue
    f_lo, f_hi = asimov_shape_interval(α, β, γ, widths, L_fb=L_fb, Aeps=Aeps, q_target=3.84)
    width = f_hi - f_lo
    label = f"W>{Wmin:.0f} GeV" if Wmax is None else f"{Wmin:.0f}<W<{Wmax:.0f} GeV"
    print(f"{label:>14}:  f ∈ [{f_lo:.3g}, {f_hi:.3g}] TeV^-4")
    if best_shape is None or width < best_shape[-1]:
        best_shape = (label, f_lo, f_hi, width)

if best_shape:
    label, f_lo, f_hi, _ = best_shape
    print(f">> Best SHAPE window: {label}  →  f ∈ [{f_lo:.3g}, {f_hi:.3g}] TeV^-4")

# =============================================================================











