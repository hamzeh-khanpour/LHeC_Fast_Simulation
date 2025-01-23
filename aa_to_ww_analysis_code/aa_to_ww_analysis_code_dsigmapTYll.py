import ROOT
import matplotlib.pyplot as plt
import numpy as np

# Matplotlib configuration for publication-quality plots
import mplhep as hep

hep.style.use("CMS")

# Function to parse the LHE file and extract lepton transverse momentum (pT) and rapidity
def parse_lhe_file(file_name):
    pt_leptons = []
    rapidities = []

    with open(file_name, "r") as file:
        in_event = False
        lepton_plus = None
        lepton_minus = None

        for line in file:
            line = line.strip()

            # Check if we are inside an <event> block
            if "<event>" in line:
                in_event = True
                lepton_plus = None
                lepton_minus = None
                continue
            if "</event>" in line:
                in_event = False
                # Calculate rapidity for the lepton pair
                if lepton_plus and lepton_minus:
                    px_pair = lepton_plus.Px() + lepton_minus.Px()
                    py_pair = lepton_plus.Py() + lepton_minus.Py()
                    pz_pair = lepton_plus.Pz() + lepton_minus.Pz()
                    energy_pair = lepton_plus.E() + lepton_minus.E()

                    if abs(energy_pair - abs(pz_pair)) > 1e-6:  # Avoid division by zero
                        rapidity = 0.5 * np.log((energy_pair + pz_pair) / (energy_pair - pz_pair))
                        rapidities.append(rapidity)
                continue

            if in_event:
                # Skip lines starting with non-numeric characters (metadata lines)
                if not line[0].isdigit() and not line[0] == "-":  # Include negative PDG IDs
                    continue

                # Split the line into components
                parts = line.split()
                if len(parts) < 10:
                    continue

                pdg_id = int(parts[0])
                px = float(parts[6])
                py = float(parts[7])
                pz = float(parts[8])
                energy = float(parts[9])

                # Check if the particle is a lepton (e+, e-, mu+, mu-, tau+, tau-)
                if pdg_id in [-11, -13, -15]:  # e+, mu+, tau+
                    lepton_plus = ROOT.TLorentzVector()
                    lepton_plus.SetPxPyPzE(px, py, pz, energy)
                elif pdg_id in [11, 13, 15]:  # e-, mu-, tau-
                    lepton_minus = ROOT.TLorentzVector()
                    lepton_minus.SetPxPyPzE(px, py, pz, energy)

                # Store transverse momentum for individual leptons
                if abs(pdg_id) in [11, 13, 15]:
                    lepton = ROOT.TLorentzVector()
                    lepton.SetPxPyPzE(px, py, pz, energy)
                    pt_leptons.append(lepton.Pt())

    return pt_leptons, rapidities

# Parse signal and background files
signal_file = "/home/hamzeh-khanpour/MG5_aMC_v3_5_7/aa_ww_semi_leptonic_NP/Events/run_01/aa_ww_semi_leptonic_NP.lhe"
background_file = "/home/hamzeh-khanpour/MG5_aMC_v3_5_7/aa_ww_semi_leptonic_SM/Events/run_01/aa_ww_semi_leptonic_SM.lhe"

pt_leptons_signal, rapidities_signal = parse_lhe_file(signal_file)
pt_leptons_background, rapidities_background = parse_lhe_file(background_file)

# Parameters for differential cross-section
signal_cross_section = 2.05983  # pb
background_cross_section = 0.0134802  # pb

num_bins = 50
pt_range = (0, 400)
rapidity_range = (-10, 10)



# Calculate bin width
bin_width_pt =  (pt_range[1] - pt_range[0]) / num_bins
bin_width_rapidity =  (rapidity_range[1] - rapidity_range[0]) / num_bins




# Normalize histograms to calculate differential cross-section
def calculate_dsigma(data, total_cross_section, num_bins, data_range):
    if len(data) == 0:
        print("Warning: No data available for differential cross-section calculation!")
        return np.zeros(num_bins), np.zeros(num_bins)  # Return empty arrays

    counts, bin_edges = np.histogram(data, bins=num_bins, range=data_range)
    bin_width = (data_range[1] - data_range[0]) / num_bins  # Calculate bin width
    dsigma = counts * (total_cross_section / len(data)) / bin_width
    return bin_edges[:-1], dsigma




# Calculate differential cross-sections
pt_bins_signal, dsigma_pt_signal = calculate_dsigma(pt_leptons_signal, signal_cross_section, num_bins, pt_range)
pt_bins_background, dsigma_pt_background = calculate_dsigma(pt_leptons_background, background_cross_section, num_bins, pt_range)

rapidity_bins_signal, dsigma_rapidity_signal = calculate_dsigma(rapidities_signal, signal_cross_section, num_bins, rapidity_range)
rapidity_bins_background, dsigma_rapidity_background = calculate_dsigma(rapidities_background, background_cross_section, num_bins, rapidity_range)



# Plot the differential cross-sections for transverse momentum
plt.figure(figsize=(10, 8))
plt.step(pt_bins_signal, dsigma_pt_signal, where="mid", alpha=0.7, label="Signal ($w^+ w^-$)", color="blue", linewidth=3)
plt.step(pt_bins_background, dsigma_pt_background, where="mid", alpha=0.7, label="Background ($w^+ w^-$)", color="red", linewidth=3)
plt.xlabel(r"$p_T^{\ell} \ [\mathrm{GeV}]$")
plt.ylabel(r"$\frac{d\sigma}{dp_T} \ [\mathrm{pb/GeV}]$")
plt.title("Differential Cross-Section: $p_T$")
plt.yscale("log")
plt.legend()
plt.grid(True, linestyle="--", alpha=0.6)
plt.tight_layout()
plt.savefig("differential_cross_section_pt.png", dpi=300)
plt.show()




# Plot the differential cross-sections for rapidity
plt.figure(figsize=(10, 8))
plt.step(rapidity_bins_signal, dsigma_rapidity_signal, where="mid", alpha=0.7, label="Signal ($w^+ w^-$)", color="blue", linewidth=3)
plt.step(rapidity_bins_background, dsigma_rapidity_background, where="mid", alpha=0.7, label="Background ($w^+ w^-$)", color="red", linewidth=3)
plt.xlabel(r"$Y^{\ell^+ \ell^-}$")
plt.ylabel(r"$\frac{d\sigma}{dY} \ [\mathrm{pb}]$")
plt.title("Differential Cross-Section: Rapidity")
plt.legend()
plt.grid(True, linestyle="--", alpha=0.6)
plt.tight_layout()
plt.savefig("differential_cross_section_rapidity.png", dpi=300)
plt.show()



