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



# Function to parse the LHE file and extract lepton transverse momentum (pT)
def parse_lhe_file(file_name):
    pt_leptons = []

    with open(file_name, "r") as file:
        in_event = False

        for line in file:
            line = line.strip()

            # Check if we are inside an <event> block
            if "<event>" in line:
                in_event = True
                continue
            if "</event>" in line:
                in_event = False
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

    return pt_leptons


    fig, ax = plt.subplots(figsize = (12.0, 10.0))
    plt.subplots_adjust(left=0.15, right=0.95, bottom=0.12, top=0.95)
    
    


# Parse signal and background files
signal_file = "/home/hamzeh-khanpour/MG5_aMC_v3_5_7/aa_ww_semi_leptonic_NP/Events/run_01/aa_ww_semi_leptonic_NP.lhe"
background_file = "/home/hamzeh-khanpour/MG5_aMC_v3_5_7/aa_ww_semi_leptonic_SM/Events/run_01/aa_ww_semi_leptonic_SM.lhe"

pt_leptons_signal = parse_lhe_file(signal_file)
pt_leptons_background = parse_lhe_file(background_file)

# Parameters for differential cross-section
signal_cross_section = 2.05983  # pb
background_cross_section = 0.0134802  # pb
num_bins = 50
pt_range = (0, 400)

# Calculate bin width
bin_width = (pt_range[1] - pt_range[0]) / num_bins

# Normalize histograms to calculate differential cross-section
def calculate_dsigma(data, total_cross_section):
    counts, bin_edges = np.histogram(data, bins=num_bins, range=pt_range)
    dsigma = counts * (total_cross_section / len(data)) / bin_width
    return bin_edges[:-1], dsigma

# Calculate differential cross-sections
pt_bins_signal, dsigma_signal = calculate_dsigma(pt_leptons_signal, signal_cross_section)
pt_bins_background, dsigma_background = calculate_dsigma(pt_leptons_background, background_cross_section)

# Plot the differential cross-sections
# plt.figure(figsize=(10, 8))

plt.step(pt_bins_signal, dsigma_signal, where="mid", alpha=0.7, label="LHeC@1.2 TeV : Signal ($w^+ w^-) [f_{M_0} / \Lambda^4$]", color="blue", linewidth=3)
plt.step(pt_bins_background, dsigma_background, where="mid", alpha=0.7, label="LHeC@1.2 TeV : SM background ($w^+ w^-$)", color="red", linewidth=3)
plt.xlabel(r"$p_T^{\ell} \ \mathrm{[GeV]}$")
plt.ylabel(r"$\frac{d\sigma}{dp_T} \ \mathrm{[pb/GeV]}$")
plt.title(r"$e^- p \to e^- w^+ w^- p \to e^- j j \ell \nu_{\ell} p$ : LHeC@1.2 TeV", fontsize=24)
plt.yscale("log")

plt.legend()

plt.grid(True, linestyle="--", alpha=0.6)
plt.tight_layout()
plt.ylim(0.000001, 1.0)


np.savetxt("dsigma_signal.txt", np.column_stack([pt_bins_signal, dsigma_signal]), header="pT [GeV], dσ/dpT [pb/GeV]")
np.savetxt("dsigma_background.txt", np.column_stack([pt_bins_background, dsigma_background]), header="pT [GeV], dσ/dpT [pb/GeV]")



plt.savefig("differential_cross_section.png", dpi=300)
plt.show()




