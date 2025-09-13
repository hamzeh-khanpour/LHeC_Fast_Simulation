
# =======================================================================
#     Hamzeh khanpour -- 2025
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
signal_file_0 = "/home/hamzeh-khanpour/MG5_aMC_v3_6_3/aa_ww_semi_leptonic_NP_1_FM0/Events/run_01/aa_ww_semi_leptonic_NP_1_FM0.lhe"
#signal_file_1 = "/home/hamzeh-khanpour/MG5_aMC_v3_5_7/aa_ww_semi_leptonic_NP_1_FM1/Events/run_01/aa_ww_semi_leptonic_NP_1_FM1.lhe"
signal_file_2 = "/home/hamzeh-khanpour/MG5_aMC_v3_6_3/aa_ww_semi_leptonic_NP_1_FM2/Events/run_02/aa_ww_semi_leptonic_NP_1_FM2_10.lhe"
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
signal_cross_section_0    = 0.0153958       # pb  FM0
#signal_cross_section_1   = 0.0       # pb  FM1
signal_cross_section_2    = 0.0153508       # pb  FM2
#signal_cross_section_3   = 0.0       # pb  FM3


background_cross_section = 0.0149219    # pb     real : 0.0154074




num_bins = 40

pt_range_lepton = (0, 300)     # Range for lepton pT
pt_range_jet = (0, 300)        # Range for leading jet pT (adjusted for higher jet momenta)
eta_range = (-10, 10)          # Range for pseudorapidity
delta_r_range = (0, 10)        # Range for Delta R is between 0 and 5
met_range = (1, 400)           # Define range for MET (adjust as needed)
centrality_range = (-10, 10)     # Centrality typically ranges from -5 to 5
exp_centrality_range = (0, 2)  # Centrality typically ranges from 0 to 2
jet_centrality_range = (0, 6) # Centrality typically ranges of jet from 0 to 10
delta_eta_jj_range  = (-10, 10)   # Pseudorapidity difference between jets from 0 to 5
m_w_hadronic_range = (1, 140)  # Range for the hadronic W boson mass
m_w_leptonic_range = (1, 140)  # Range for the leptonic W boson mass

# >>> ADD THESE LINES FOR RECO M_WW <<<
m_w_hadronic_leptonic_range = (165, 1000)  # WW invariant mass; threshold ~ 2*mW up to ~LHeC √s_γγ





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






# Normalize histograms to calculate differential cross-section (width-aware)
def calculate_dsigma(data, total_cross_section, bin_width, data_range):
    """
    Returns (x_array, dsigma) where x_array are the *left bin edges*.
    NOTE: We ignore the scalar `bin_width` you pass in and instead divide by the
    true per-bin widths from the histogram edges. This stays correct if `num_bins`
    changes or if you ever switch to non-uniform bins.

    If you plot with plt.step(..., where="mid"), use bin *centers* instead:
        edges = x_array
        widths = np.diff(np.r_[edges, edges[-1] + (edges[-1]-edges[-2])])
        centers = edges + 0.5*widths
        plt.step(centers, dsigma, where="mid", ...)
    Otherwise, keep left edges and use where="post".
    """
    counts, edges = np.histogram(data, bins=num_bins, range=data_range)

    # Guard against empty selections
    n_total = max(len(data), 1)

    # Per-bin widths (handles any binning)
    widths = np.diff(edges).astype(float)

    # dσ/dx_i = counts_i * (σ_tot / N_events) / Δx_i
    with np.errstate(divide="ignore", invalid="ignore"):
        dsigma = counts * (total_cross_section / n_total) / widths
        # Replace any NaN/inf (can only happen with zero-width bins) with 0
        dsigma = np.nan_to_num(dsigma, nan=0.0, posinf=0.0, neginf=0.0)

    # Return left edges to preserve your current calling convention
    return edges[:-1], dsigma






# >>> ADD THESE LINES FOR RECO M_WW <<<
m_w_hadronic_leptonic_bins_signal_0, dsigma_signal_m_w_hadronic_leptonic_0 = calculate_dsigma(    m_w_hadronic_leptonic_signal_0, signal_cross_section_0, bin_width_m_w_hadronic_leptonic, m_w_hadronic_leptonic_range)
m_w_hadronic_leptonic_bins_signal_2, dsigma_signal_m_w_hadronic_leptonic_2 = calculate_dsigma(    m_w_hadronic_leptonic_signal_2, signal_cross_section_2, bin_width_m_w_hadronic_leptonic, m_w_hadronic_leptonic_range)
m_w_hadronic_leptonic_bins_background, dsigma_background_m_w_hadronic_leptonic = calculate_dsigma(    m_w_hadronic_leptonic_background, background_cross_section, bin_width_m_w_hadronic_leptonic, m_w_hadronic_leptonic_range)








# =======================================================================




# --- Normalize the reconstructed M_WW distributions ---
s0 = np.sum(dsigma_signal_m_w_hadronic_leptonic_0)
s2 = np.sum(dsigma_signal_m_w_hadronic_leptonic_2)
sb = np.sum(dsigma_background_m_w_hadronic_leptonic)

dsigma_signal_m_w_hadronic_leptonic_0_norm = dsigma_signal_m_w_hadronic_leptonic_0 / s0 if s0 > 0 else dsigma_signal_m_w_hadronic_leptonic_0
dsigma_signal_m_w_hadronic_leptonic_2_norm = dsigma_signal_m_w_hadronic_leptonic_2 / s2 if s2 > 0 else dsigma_signal_m_w_hadronic_leptonic_2
dsigma_background_m_w_hadronic_leptonic_norm = dsigma_background_m_w_hadronic_leptonic / sb if sb > 0 else dsigma_background_m_w_hadronic_leptonic




# --------------------------------------------
# Plot the normalized Reconstructed WW Invariant Mass (M_WW)
plt.step(m_w_hadronic_leptonic_bins_signal_0, dsigma_signal_m_w_hadronic_leptonic_0_norm, where="mid", alpha=0.7, label="LHeC@1.2 TeV : Signal ($w^+ w^-) [f_{M_0} / \\Lambda^4$)", color="red", linewidth=3)
plt.step(m_w_hadronic_leptonic_bins_signal_2, dsigma_signal_m_w_hadronic_leptonic_2_norm, where="mid", alpha=0.7, label="LHeC@1.2 TeV : Signal ($w^+ w^-) [f_{M_2} / \\Lambda^4$)", color="green", linewidth=3)
plt.step(m_w_hadronic_leptonic_bins_background, dsigma_background_m_w_hadronic_leptonic_norm, where="mid", alpha=0.7, label="LHeC@1.2 TeV : SM background ($w^+ w^-$)", color="blue", linewidth=3)
plt.xlabel(r"$M_{WW} \ \mathrm{[GeV]}$")
plt.ylabel("Normalized Distribution")
plt.title(r"$e^- p \to e^- w^+ w^- p \to e^- j j \ell \nu_{\ell} p$ : LHeC@1.2 TeV", fontsize=24)
# plt.yscale("log")  # usually off for normalized shapes; keep commented like your leptonic plot
plt.legend()
plt.grid(True, linestyle="--", alpha=0.6)
plt.tight_layout()
plt.savefig("normalized_cross_section_m_ww_reco.pdf", dpi=600)
plt.show()




# --------------------------------------------
# Plot the differential cross-sections for Reconstructed WW Invariant Mass (M_WW)
plt.step(m_w_hadronic_leptonic_bins_signal_0, dsigma_signal_m_w_hadronic_leptonic_0, where="mid", alpha=0.7, label="LHeC@1.2 TeV : Signal ($w^+ w^-) [f_{M_0} / \\Lambda^4$)", color="red", linewidth=3)
plt.step(m_w_hadronic_leptonic_bins_signal_2, dsigma_signal_m_w_hadronic_leptonic_2, where="mid", alpha=0.7, label="LHeC@1.2 TeV : Signal ($w^+ w^-) [f_{M_2} / \\Lambda^4$)", color="green", linewidth=3)
plt.step(m_w_hadronic_leptonic_bins_background, dsigma_background_m_w_hadronic_leptonic, where="mid", alpha=0.7, label="LHeC@1.2 TeV : SM background ($w^+ w^-$)", color="blue", linewidth=3)
plt.xlabel(r"$M_{WW} \ \mathrm{[GeV]}$")
plt.ylabel(r"$\frac{d\sigma}{dM_{WW}} \ \mathrm{[pb/GeV]}$")
plt.title(r"$e^- p \to e^- w^+ w^- p \to e^- j j \ell \nu_{\ell} p$ : LHeC@1.2 TeV", fontsize=24)
plt.yscale("log")
plt.ylim(0.000001, 0.01)
plt.legend()
plt.grid(True, linestyle="--", alpha=0.6)
plt.tight_layout()
plt.savefig("differential_cross_section_m_ww_reco.pdf", dpi=600)
plt.show()









# =======================================================================






## Save data for pT lepton
#np.savetxt("dsigma_signal_pt_0.txt", np.column_stack([pt_bins_signal_0, dsigma_signal_pt_0]), header="pT [GeV], dσ/dpT [pb/GeV]")
#np.savetxt("dsigma_signal_pt_2.txt", np.column_stack([pt_bins_signal_2, dsigma_signal_pt_2]), header="pT [GeV], dσ/dpT [pb/GeV]")
#np.savetxt("dsigma_background_pt.txt", np.column_stack([pt_bins_background, dsigma_background_pt]), header="pT [GeV], dσ/dpT [pb/GeV]")


## Save data for eta lepton
#np.savetxt("dsigma_signal_eta_0.txt", np.column_stack([eta_bins_signal_0, dsigma_signal_eta_0]), header="eta, dσ/deta [pb]")
#np.savetxt("dsigma_signal_eta_2.txt", np.column_stack([eta_bins_signal_2, dsigma_signal_eta_2]), header="eta, dσ/deta [pb]")
#np.savetxt("dsigma_background_eta.txt", np.column_stack([eta_bins_background, dsigma_background_eta]), header="eta, dσ/deta [pb]")


## Save data for leading jet pT
#np.savetxt("dsigma_signal_leading_jet_pt_0.txt", np.column_stack([pt_jet_bins_signal_0, dsigma_signal_jet_pt_0]), header="pT [GeV], dσ/dpT [pb/GeV]")
#np.savetxt("dsigma_signal_leading_jet_pt_2.txt", np.column_stack([pt_jet_bins_signal_2, dsigma_signal_jet_pt_2]), header="pT [GeV], dσ/dpT [pb/GeV]")
#np.savetxt("dsigma_background_leading_jet_pt.txt", np.column_stack([pt_jet_bins_background, dsigma_background_jet_pt]), header="pT [GeV], dσ/dpT [pb/GeV]")


## Save data for Delta R
#np.savetxt("dsigma_signal_delta_r_0.txt", np.column_stack([delta_r_bins_signal_0, dsigma_signal_delta_r_0]), header="Delta R, dσ/d(Delta R) [pb]")
#np.savetxt("dsigma_signal_delta_r_2.txt", np.column_stack([delta_r_bins_signal_2, dsigma_signal_delta_r_2]), header="Delta R, dσ/d(Delta R) [pb]")
#np.savetxt("dsigma_background_delta_r.txt", np.column_stack([delta_r_bins_background, dsigma_background_delta_r]), header="Delta R, dσ/d(Delta R) [pb]")


## Save data for Missing Transverse Energy (MET)
#np.savetxt("dsigma_signal_met_0.txt", np.column_stack([met_bins_signal_0, dsigma_signal_met_0]), header="MET [GeV], dσ/d(MET) [pb/GeV]")
#np.savetxt("dsigma_signal_met_2.txt", np.column_stack([met_bins_signal_2, dsigma_signal_met_2]), header="MET [GeV], dσ/d(MET) [pb/GeV]")
#np.savetxt("dsigma_background_met.txt", np.column_stack([met_bins_background, dsigma_background_met]), header="MET [GeV], dσ/d(MET) [pb/GeV]")


## Save data for Lepton Centrality
#np.savetxt("dsigma_signal_centrality_0.txt", np.column_stack([centrality_bins_signal_0, dsigma_signal_centrality_0]), header="C_l, dσ/d(C_l) [pb]")
#np.savetxt("dsigma_signal_centrality_2.txt", np.column_stack([centrality_bins_signal_2, dsigma_signal_centrality_2]), header="C_l, dσ/d(C_l) [pb]")
#np.savetxt("dsigma_background_centrality.txt", np.column_stack([centrality_bins_background, dsigma_background_centrality]), header="C_l, dσ/d(C_l) [pb]")


## Save data for Exponential Lepton Centrality (exp_centrality)
#np.savetxt("dsigma_signal_exp_centrality_0.txt", np.column_stack([exp_centrality_bins_signal_0, dsigma_signal_exp_centrality_0]), header="C_l_exp, dσ/d(C_l_exp) [pb]")
#np.savetxt("dsigma_signal_exp_centrality_2.txt", np.column_stack([exp_centrality_bins_signal_2, dsigma_signal_exp_centrality_2]), header="C_l_exp, dσ/d(C_l_exp) [pb]")
#np.savetxt("dsigma_background_exp_centrality.txt", np.column_stack([exp_centrality_bins_background, dsigma_background_exp_centrality]), header="C_l_exp, dσ/d(C_l_exp) [pb]")


## Save data for Jet Centrality
#np.savetxt("dsigma_signal_jet_centrality_0.txt", np.column_stack([jet_centrality_bins_signal_0, dsigma_signal_jet_centrality_0]), header="C_j, dσ/d(C_j) [pb]")
#np.savetxt("dsigma_signal_jet_centrality_2.txt", np.column_stack([jet_centrality_bins_signal_2, dsigma_signal_jet_centrality_2]), header="C_j, dσ/d(C_j) [pb]")
#np.savetxt("dsigma_background_jet_centrality.txt", np.column_stack([jet_centrality_bins_background, dsigma_background_jet_centrality]), header="C_j, dσ/d(C_j) [pb]")



## Save data for Pseudorapidity Difference Between Jets (Δηjj)
#np.savetxt("dsigma_signal_delta_eta_jj_0.txt", np.column_stack([delta_eta_jj_bins_signal_0, dsigma_signal_delta_eta_jj_0]), header="Δη_jj, dσ/d(Δη_jj) [pb]")
#np.savetxt("dsigma_signal_delta_eta_jj_2.txt", np.column_stack([delta_eta_jj_bins_signal_2, dsigma_signal_delta_eta_jj_2]), header="Δη_jj, dσ/d(Δη_jj) [pb]")
#np.savetxt("dsigma_background_delta_eta_jj.txt", np.column_stack([delta_eta_jj_bins_background, dsigma_background_delta_eta_jj]), header="Δη_jj, dσ/d(Δη_jj) [pb]")






## Save data for Hadronic W Boson Invariant Mass Distribution (M_W^{jj})
#np.savetxt("dsigma_signal_m_w_hadronic_0.txt",   np.column_stack([m_w_hadronic_bins_signal_0, dsigma_signal_m_w_hadronic_0]), header="M_W^{jj} [GeV]   dσ/dM_W^{jj} [pb/GeV]", fmt="%.6e")
#np.savetxt("dsigma_signal_m_w_hadronic_2.txt",   np.column_stack([m_w_hadronic_bins_signal_2, dsigma_signal_m_w_hadronic_2]), header="M_W^{jj} [GeV]   dσ/dM_W^{jj} [pb/GeV]", fmt="%.6e")
#np.savetxt("dsigma_background_m_w_hadronic.txt", np.column_stack([m_w_hadronic_bins_background, dsigma_background_m_w_hadronic]), header="M_W^{jj} [GeV]   dσ/dM_W^{jj} [pb/GeV]", fmt="%.6e")



## Save data for Hadronic W Boson Invariant Mass Distribution (M_W^{jj})
#np.savetxt("dsigma_signal_m_w_leptonic_0.txt",   np.column_stack([m_w_leptonic_bins_signal_0, dsigma_signal_m_w_leptonic_0]), header="M_W^{jj} [GeV]   dσ/dM_W^{jj} [pb/GeV]", fmt="%.6e")
#np.savetxt("dsigma_signal_m_w_leptonic_2.txt",   np.column_stack([m_w_leptonic_bins_signal_2, dsigma_signal_m_w_leptonic_2]), header="M_W^{jj} [GeV]   dσ/dM_W^{jj} [pb/GeV]", fmt="%.6e")
#np.savetxt("dsigma_background_m_w_leptonic.txt", np.column_stack([m_w_leptonic_bins_background, dsigma_background_m_w_leptonic]), header="M_W^{jj} [GeV]   dσ/dM_W^{jj} [pb/GeV]", fmt="%.6e")




