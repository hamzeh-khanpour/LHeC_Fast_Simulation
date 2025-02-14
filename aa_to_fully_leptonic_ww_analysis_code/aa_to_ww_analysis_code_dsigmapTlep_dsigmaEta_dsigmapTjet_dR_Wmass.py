
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


# Function to parse the LHE file and extract lepton transverse momentum (pT), eta, and Delta R
# Jets are vetoed, meaning only pure leptonic events are considered.

def parse_lhe_file(file_name):
    pt_leptons = []
    eta_leptons = []
    delta_r_values = []  # Delta R between leptons
    missing_transverse_energy = []  # MET (Missing Transverse Energy)
    rapidity_lepton_pair = []  # Rapidity of lepton pair

    with open(file_name, "r") as file:
        in_event = False
        lepton_plus = None
        lepton_minus = None

        for line in file:
            line = line.strip()

            # Check if we are inside an <event> block
            if "<event>" in line:
                in_event = True
                leptons = []  # Store lepton 4-momenta
                neutrinos = []  # Store neutrino 4-momenta for MET calculation
                jets = []  # Keep track of jets for vetoing
                lepton_plus = None
                lepton_minus = None  # Reset lepton tracking
                continue
            if "</event>" in line:
                in_event = False

                # Veto events that contain any jets
                if len(jets) > 0:
                    continue  # Skip this event since it contains jets

                # Compute MET from neutrinos' pT
                met_x = sum(nu.Px() for nu in neutrinos)
                met_y = sum(nu.Py() for nu in neutrinos)
                met = np.sqrt(met_x**2 + met_y**2)
                missing_transverse_energy.append(met)

                # Compute ΔR between leptons
                if len(leptons) >= 2:
                    for i in range(len(leptons)):
                        for j in range(i + 1, len(leptons)):
                            delta_r = leptons[i].DeltaR(leptons[j])
                            delta_r_values.append(delta_r)

                # Compute rapidity of the lepton pair if both leptons exist
                if lepton_plus and lepton_minus:
                    px_pair = lepton_plus["px"] + lepton_minus["px"]
                    py_pair = lepton_plus["py"] + lepton_minus["py"]
                    pz_pair = lepton_plus["pz"] + lepton_minus["pz"]
                    energy_pair = lepton_plus["energy"] + lepton_minus["energy"]

                    # Compute rapidity and store it
                    if abs(energy_pair - abs(pz_pair)) > 1e-6:
                        rapidity = 0.5 * np.log((energy_pair + pz_pair) / (energy_pair - pz_pair))
                        rapidity_lepton_pair.append(rapidity)

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
                if abs(pdg_id) in [11, 13, 15]:  # e, μ, τ
                    lepton = ROOT.TLorentzVector()
                    lepton.SetPxPyPzE(px, py, pz, energy)

                    # Store lepton properties
                    leptons.append(lepton)
                    pt_leptons.append(lepton.Pt())
                    eta_leptons.append(lepton.Eta())

                    # Identify lepton charge and store information
                    lepton_info = {"px": px, "py": py, "pz": pz, "energy": energy}

                    if pdg_id > 0:
                        lepton_plus = lepton_info
                    elif pdg_id < 0:
                        lepton_minus = lepton_info

                # Check if the particle is a neutrino (νe, νμ, ντ)
                if abs(pdg_id) in [12, 14, 16]:  # Neutrinos
                    neutrino = ROOT.TLorentzVector()
                    neutrino.SetPxPyPzE(px, py, pz, energy)

                    # Store neutrino 4-momentum for MET calculation
                    neutrinos.append(neutrino)

                # Check if the particle is a jet (not a lepton, neutrino, or photon)
                if abs(pdg_id) not in [11, 12, 13, 14, 15, 16, 22] and int(parts[1]) == 1:
                    jet = ROOT.TLorentzVector()
                    jet.SetPxPyPzE(px, py, pz, energy)

                    # Store jet information (for vetoing)
                    jets.append(jet)

    return pt_leptons, eta_leptons, delta_r_values, missing_transverse_energy, rapidity_lepton_pair




# =======================================================================




# Adjust subplot layout
plt.subplots_adjust(left=0.15, right=0.95, bottom=0.12, top=0.95)



# Define LHE file paths for signals and background
signal_file_0 = "/home/hamzeh-khanpour/MG5_aMC_v3_5_7/aa_ww_fully_leptonic_NP_FM0/Events/run_01/aa_ww_fully_leptonic_NP_FM0.lhe"
signal_file_2 = "/home/hamzeh-khanpour/MG5_aMC_v3_5_7/aa_ww_fully_leptonic_NP_FM2/Events/run_01/aa_ww_fully_leptonic_NP_FM2.lhe"
background_file = "/home/hamzeh-khanpour/MG5_aMC_v3_5_7/aa_ww_fully_leptonic_SM/Events/run_01/aa_ww_fully_leptonic_SM.lhe"





# Parse LHE files, keeping only leptonic events (jets are vetoed)
(pt_leptons_signal_0, eta_leptons_signal_0, delta_r_signal_0, met_signal_0, rapidity_lepton_pair_signal_0) = parse_lhe_file(signal_file_0)
(pt_leptons_signal_2, eta_leptons_signal_2, delta_r_signal_2, met_signal_2, rapidity_lepton_pair_signal_2) = parse_lhe_file(signal_file_2)
(pt_leptons_background, eta_leptons_background, delta_r_background, met_background, rapidity_lepton_pair_background) = parse_lhe_file(background_file)



# =======================================================================



# Parameters for differential cross-section

signal_cross_section_0    =  0.503747      # pb  FM0
#signal_cross_section_1   =  0.0451322    # pb  FM1
signal_cross_section_2    =  19.66510      # pb  FM2
#signal_cross_section_3   =  1.505529     # pb  FM3


background_cross_section  =  0.00351599    # pb


num_bins = 50

pt_range_lepton = (0, 500)     # Range for lepton pT
eta_range = (-10, 10)          # Range for pseudorapidity
delta_r_range = (0, 10)        # Range for Delta R between leptons
met_range = (1, 500)           # Define range for MET (adjust as needed)
rapidity_range = (-5, 5)       # Define range for lepton pair rapidity



# Calculate bin widths
bin_width_pt_lepton = (pt_range_lepton[1] - pt_range_lepton[0]) / num_bins
bin_width_eta = (eta_range[1] - eta_range[0]) / num_bins
bin_width_delta_r = (delta_r_range[1] - delta_r_range[0]) / num_bins
bin_width_met = (met_range[1] - met_range[0]) / num_bins  # Bin width for MET
bin_width_rapidity = (rapidity_range[1] - rapidity_range[0]) / num_bins  # Bin width for lepton pair rapidity





# Normalize histograms to calculate differential cross-section
def calculate_dsigma(data, total_cross_section, bin_width, data_range):
    counts, bin_edges = np.histogram(data, bins=num_bins, range=data_range)
    dsigma = counts * (total_cross_section / len(data)) / bin_width
    return bin_edges[:-1], dsigma



# Calculate differential cross-sections for lepton pT and eta
pt_bins_signal_0, dsigma_signal_pt_0 = calculate_dsigma(pt_leptons_signal_0, signal_cross_section_0, bin_width_pt_lepton, pt_range_lepton)
pt_bins_signal_2, dsigma_signal_pt_2 = calculate_dsigma(pt_leptons_signal_2, signal_cross_section_2, bin_width_pt_lepton, pt_range_lepton)
pt_bins_background, dsigma_background_pt = calculate_dsigma(pt_leptons_background, background_cross_section, bin_width_pt_lepton, pt_range_lepton)

eta_bins_signal_0, dsigma_signal_eta_0 = calculate_dsigma(eta_leptons_signal_0, signal_cross_section_0, bin_width_eta, eta_range)
eta_bins_signal_2, dsigma_signal_eta_2 = calculate_dsigma(eta_leptons_signal_2, signal_cross_section_2, bin_width_eta, eta_range)
eta_bins_background, dsigma_background_eta = calculate_dsigma(eta_leptons_background, background_cross_section, bin_width_eta, eta_range)

# Normalize Delta R histograms
delta_r_bins_signal_0, dsigma_signal_delta_r_0 = calculate_dsigma(delta_r_signal_0, signal_cross_section_0, bin_width_delta_r, delta_r_range)
delta_r_bins_signal_2, dsigma_signal_delta_r_2 = calculate_dsigma(delta_r_signal_2, signal_cross_section_2, bin_width_delta_r, delta_r_range)
delta_r_bins_background, dsigma_background_delta_r = calculate_dsigma(delta_r_background, background_cross_section, bin_width_delta_r, delta_r_range)


# Calculate differential cross-sections for MET
met_bins_signal_0, dsigma_signal_met_0 = calculate_dsigma(met_signal_0, signal_cross_section_0, bin_width_met, met_range)
met_bins_signal_2, dsigma_signal_met_2 = calculate_dsigma(met_signal_2, signal_cross_section_2, bin_width_met, met_range)
met_bins_background, dsigma_background_met = calculate_dsigma(met_background, background_cross_section, bin_width_met, met_range)



# Calculate differential cross-sections for Lepton Pair Rapidity
rapidity_bins_signal_0, dsigma_signal_rapidity_0 = calculate_dsigma(rapidity_lepton_pair_signal_0, signal_cross_section_0, bin_width_rapidity, rapidity_range)
rapidity_bins_signal_2, dsigma_signal_rapidity_2 = calculate_dsigma(rapidity_lepton_pair_signal_2, signal_cross_section_2, bin_width_rapidity, rapidity_range)
rapidity_bins_background, dsigma_background_rapidity = calculate_dsigma(rapidity_lepton_pair_background, background_cross_section, bin_width_rapidity, rapidity_range)




# =======================================================================




# Plot the differential cross-sections
# plt.figure(figsize=(10, 8))

plt.step(pt_bins_signal_0, dsigma_signal_pt_0, where="mid", alpha=0.7, label="LHeC@1.2 TeV : Signal ($w^+ w^-) [f_{M_0} / \Lambda^4$]", color="red", linewidth=3)
plt.step(pt_bins_signal_2, dsigma_signal_pt_2, where="mid", alpha=0.7, label="LHeC@1.2 TeV : Signal ($w^+ w^-) [f_{M_2} / \Lambda^4$]", color="green", linewidth=3)
plt.step(pt_bins_background, dsigma_background_pt, where="mid", alpha=0.7, label="LHeC@1.2 TeV : SM background ($w^+ w^-$)", color="blue", linewidth=3)
plt.xlabel(r"$p_T^{\ell} \ \mathrm{[GeV]}$")
plt.ylabel(r"$\frac{d\sigma}{dp_T^{\ell}} \ \mathrm{[pb/GeV]}$")
plt.title(r"$e^- p \to e^- W^+ W^- p \to e^- \ell^+ \nu_{\ell} \ell^- \bar{\nu}_{\ell} p$ : LHeC@1.2 TeV", fontsize=24)
plt.yscale("log")
plt.legend()
plt.grid(True, linestyle="--", alpha=0.6)
plt.tight_layout()
plt.ylim(0.00001, 10.0)
plt.savefig("differential_cross_section_pt_fully_leptonic.png", dpi=600)
plt.show()




# Plot the differential cross-sections for eta
#plt.figure(figsize=(10, 8))

plt.step(eta_bins_signal_0, dsigma_signal_eta_0, where="mid", alpha=0.7, label="LHeC@1.2 TeV : Signal ($w^+ w^-) [f_{M_0} / \Lambda^4$]", color="red", linewidth=3)
plt.step(eta_bins_signal_2, dsigma_signal_eta_2, where="mid", alpha=0.7, label="LHeC@1.2 TeV : Signal ($w^+ w^-) [f_{M_2} / \Lambda^4$]", color="green", linewidth=3)
plt.step(eta_bins_background, dsigma_background_eta, where="mid", alpha=0.7, label="LHeC@1.2 TeV : SM background ($w^+ w^-$)", color="blue", linewidth=3)
plt.xlabel(r"$\eta^{\ell}$")
plt.ylabel(r"$\frac{d\sigma}{d\eta^{\ell}} \ \mathrm{[pb]}$")
plt.title(r"$e^- p \to e^- W^+ W^- p \to e^- \ell^+ \nu_{\ell} \ell^- \bar{\nu}_{\ell} p$ : LHeC@1.2 TeV", fontsize=24)
plt.yscale("log")
plt.ylim(0.0001, 100.0)
plt.legend()
plt.grid(True, linestyle="--", alpha=0.6)
plt.tight_layout()
plt.savefig("differential_cross_section_eta_fully_leptonic.png", dpi=600)
plt.show()





# Plot Delta R (between leptons) differential cross-section
plt.step(delta_r_bins_signal_0, dsigma_signal_delta_r_0, where="mid", alpha=0.7, label=r"LHeC@1.2 TeV : Signal ($W^+ W^-$) [$f_{M_0} / \Lambda^4$]", color="red", linewidth=3)
plt.step(delta_r_bins_signal_2, dsigma_signal_delta_r_2, where="mid", alpha=0.7, label=r"LHeC@1.2 TeV : Signal ($W^+ W^-$) [$f_{M_2} / \Lambda^4$]", color="green", linewidth=3)
plt.step(delta_r_bins_background, dsigma_background_delta_r, where="mid", alpha=0.7, label=r"LHeC@1.2 TeV : SM background ($W^+ W^-$)", color="blue", linewidth=3)
# Updated X-axis label for correct interpretation
plt.xlabel(r"$\Delta R(\ell_1, \ell_2)$")
# Y-axis label remains the same
plt.ylabel(r"$\frac{d\sigma}{d\Delta R} \ \mathrm{[pb]}$")
# Updated plot title to reflect the fully leptonic process
plt.title(r"$e^- p \to e^- W^+ W^- p \to e^- \ell^+ \nu_{\ell} \ell^- \bar{\nu}_{\ell} p$ : LHeC@1.2 TeV", fontsize=24)
plt.yscale("log")
plt.ylim(0.0001, 1000.0)
plt.legend()
plt.grid(True, linestyle="--", alpha=0.6)
plt.tight_layout()
plt.savefig("differential_cross_section_delta_r_fully_leptonic.png", dpi=600)
plt.show()





# Plot the differential cross-sections for Missing Transverse Energy (MET)
plt.step(met_bins_signal_0, dsigma_signal_met_0, where="mid", alpha=0.7, label="LHeC@1.2 TeV : Signal ($w^+ w^-) [f_{M_0} / \Lambda^4$]", color="red", linewidth=3)
plt.step(met_bins_signal_2, dsigma_signal_met_2, where="mid", alpha=0.7, label="LHeC@1.2 TeV : Signal ($w^+ w^-) [f_{M_2} / \Lambda^4$]", color="green", linewidth=3)
plt.step(met_bins_background, dsigma_background_met, where="mid", alpha=0.7, label="LHeC@1.2 TeV : SM background ($w^+ w^-$)", color="blue", linewidth=3)
plt.xlabel(r"$\mathrm{MET} \ \mathrm{[GeV]}$")
plt.ylabel(r"$\frac{d\sigma}{d\mathrm{MET}} \ \mathrm{[pb/GeV]}$")
plt.title(r"$e^- p \to e^- W^+ W^- p \to e^- \ell^+ \nu_{\ell} \ell^- \bar{\nu}_{\ell} p$ : LHeC@1.2 TeV", fontsize=24)
plt.yscale("log")
plt.legend()
plt.grid(True, linestyle="--", alpha=0.6)
plt.tight_layout()
plt.ylim(0.000001, 1.0)
plt.savefig("differential_cross_section_met_fully_leptonic.png", dpi=600)
plt.show()






# Plot the differential cross-sections for Rapidity of Lepton Pair
plt.step(rapidity_bins_signal_0, dsigma_signal_rapidity_0, where="mid", alpha=0.7, label="LHeC@1.2 TeV : Signal ($W^+ W^-$) [$f_{M_0} / \Lambda^4$]", color="red", linewidth=3)
plt.step(rapidity_bins_signal_2, dsigma_signal_rapidity_2, where="mid", alpha=0.7, label="LHeC@1.2 TeV : Signal ($W^+ W^-$) [$f_{M_2} / \Lambda^4$]", color="green", linewidth=3)
plt.step(rapidity_bins_background, dsigma_background_rapidity, where="mid", alpha=0.7, label="LHeC@1.2 TeV : SM background ($W^+ W^-$)", color="blue", linewidth=3)
plt.xlabel(r"$Y(\ell^+ \ell^-)$")
plt.ylabel(r"$\frac{d\sigma}{dY(\ell^+ \ell^-)} \ \mathrm{[pb]}$")
plt.title(r"$e^- p \to e^- W^+ W^- p \to e^- \ell^+ \nu_{\ell} \ell^- \bar{\nu}_{\ell} p$ : LHeC@1.2 TeV", fontsize=24)
plt.yscale("log")
plt.legend()
plt.grid(True, linestyle="--", alpha=0.6)
plt.tight_layout()
plt.ylim(0.0001, 100.0)
plt.savefig("differential_cross_section_rapidity_fully_leptonic.png", dpi=600)
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



