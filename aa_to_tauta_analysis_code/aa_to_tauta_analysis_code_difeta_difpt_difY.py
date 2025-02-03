

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


# Function to parse the LHE file and extract Pt, eta of tau+, tau-, and rapidity of tau+ tau- pair
def parse_lhe_file(file_name):
    pt_tau_plus = []
    eta_tau_plus = []
    pt_tau_minus = []
    eta_tau_minus = []
    rapidity_tau_pair = []  # To store rapidity of tau+ tau- pair

    with open(file_name, "r") as file:
        in_event = False
        tau_plus = None
        tau_minus = None

        for line in file:
            line = line.strip()
            # Check if we are inside an <event> block
            if "<event>" in line:
                in_event = True
                tau_plus = None
                tau_minus = None
                continue
            if "</event>" in line:
                in_event = False

                # Calculate rapidity of the tau+ tau- pair if both are present
                if tau_plus and tau_minus:
                    px_pair = tau_plus["px"] + tau_minus["px"]
                    py_pair = tau_plus["py"] + tau_minus["py"]
                    pz_pair = tau_plus["pz"] + tau_minus["pz"]
                    energy_pair = tau_plus["energy"] + tau_minus["energy"]

                    if abs(energy_pair - abs(pz_pair)) > 1e-6:  # Avoid division by zero
                        rapidity = 0.5 * np.log((energy_pair + pz_pair) / (energy_pair - pz_pair))
                        rapidity_tau_pair.append(rapidity)

                continue

            if in_event:
                # Skip lines starting with non-numeric characters (metadata lines)
                if not line[0].isdigit() and not line[0] == '-':  # Include negative PDG IDs
                    continue

                # Split the line into components
                parts = line.split()
                # Check if the line has enough elements to extract particle data
                if len(parts) < 10:
                    continue

                try:
                    pdg_id = int(parts[0])
                    px = float(parts[6])
                    py = float(parts[7])
                    pz = float(parts[8])
                    energy = float(parts[9])

                    # Calculate transverse momentum Pt = sqrt(px^2 + py^2)
                    pt = (px**2 + py**2)**0.5

                    # Calculate pseudorapidity eta = 0.5 * ln((E + pz) / (E - pz)), avoid division by zero
                    if abs(energy - abs(pz)) > 1e-6:  # Add small epsilon to avoid precision issues
                        eta = 0.5 * np.log((energy + pz) / (energy - pz))

                        # Check if the particle is tau+ (PDG ID = 15)
                        if pdg_id == 15:
                            pt_tau_plus.append(pt)
                            eta_tau_plus.append(eta)
                            tau_plus = {"px": px, "py": py, "pz": pz, "energy": energy}

                        # Check if the particle is tau- (PDG ID = -15)
                        if pdg_id == -15:
                            pt_tau_minus.append(pt)
                            eta_tau_minus.append(eta)
                            tau_minus = {"px": px, "py": py, "pz": pz, "energy": energy}

                except ValueError as e:
                    print(f"Error parsing line: {line}, error: {e}")
                    continue

    return pt_tau_plus, eta_tau_plus, pt_tau_minus, eta_tau_minus, rapidity_tau_pair

    fig, ax = plt.subplots(figsize = (8.0, 9.0))
    plt.subplots_adjust(left=0.15, right=0.95, bottom=0.12, top=0.95)


# Function to calculate and plot the rapidity distribution with weights
def plot_weighted_distribution(data, bins, range, color, xlabel, ylabel, title, filename, label, integrated_cross_section, integrated_luminosity):
    num_entries = len(data)
    bin_width = (range[1] - range[0]) / bins  # Calculate bin width
    event_weight = (integrated_cross_section * integrated_luminosity) / num_entries  # Weight per event in pb

#    plt.figure(figsize=(8, 9))

    counts, bin_edges, _ = plt.hist(data, bins=bins, range=range, color=color, alpha=0.7,
                                    label=f"MadGraph5_aMC@NLO", weights=[event_weight] * num_entries, edgecolor="black")
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.legend()
    plt.grid(True, linestyle="--")
    plt.tight_layout()
    plt.ylim(0.0001, 12.0)



    ## Add annotations for bin values
    #bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    #for x, y in zip(bin_centers, counts):
        #if y > 0:
            #plt.text(x, y, f"{y:.2f}", ha="center")

    plt.savefig(filename, dpi=300)
    plt.show()


# Parse the file to get Pt, eta of tau+ and tau-, and rapidity of tau+ tau- pair
file_name = "/home/hamzeh-khanpour/MG5_aMC_v3_5_7/aa_tautau_SM/Events/run_01/aa_to_tautau_SM.lhe"
pt_tau_plus, eta_tau_plus, pt_tau_minus, eta_tau_minus, rapidity_tau_pair = parse_lhe_file(file_name)

# Define integrated luminosity and cross-section
integrated_luminosity = 1.0  # fb^-1
integrated_cross_section = 46.841  # pb

# Plot the weighted rapidity distribution for tau+ tau- pair
plot_weighted_distribution(rapidity_tau_pair, bins=20, range=(-10, 10), color="purple",
                           xlabel=r"$Y_{\tau^+ \tau^-}$", ylabel=r"$\frac{d\sigma}{dY^{\tau^+ \tau^-}}$ [pb/GeV]",
                           title=r"$e^- p \to e^- \tau^+ \tau^- p$ (LHeC@1.2TeV) ", filename="rapidity_tau_pair_weighted_distribution.png",
                           label=r"$Y_{\tau^+ \tau^-}$", integrated_cross_section=integrated_cross_section,
                           integrated_luminosity=integrated_luminosity) 



